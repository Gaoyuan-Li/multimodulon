"""
Multi-view ICA implementation for MultiModulon analysis.

This module provides GPU-accelerated multi-view ICA functionality
using PyTorch natively.
"""

import sys
import time
from pathlib import Path
from typing import Dict, Optional, Tuple, List

import numpy as np
import pandas as pd

from .multiview_ica_optimization import calculate_cohens_d_effect_size, calculate_average_effect_sizes

try:
    import torch
    import torch.nn as nn
    import geotorch
    import scipy.optimize
    from torch.utils.data import DataLoader, TensorDataset
except ImportError as e:
    raise ImportError(
        "PyTorch and geotorch are required for multi-view ICA. "
        "Please install them with: pip install torch==2.6.0 torchvision==0.21.0 "
        "torchaudio==2.6.0 geotorch==0.3.0 --index-url https://download.pytorch.org/whl/cu124"
    ) from e


class MSIICA(nn.Module):
    """Multi-view ICA model using orthogonal constraints."""
    
    def __init__(self, n_in, n_out, U=None, ortho=True):
        super().__init__()
        self.W = nn.Linear(n_in, n_out, bias=False)
        geotorch.orthogonal(self.W, "weight")
        if U is not None:
            self.W.weight = U.contiguous()
                
    def forward(self, Xw):
        S = self.W(Xw)
        return S


def whiten(X, rank):
    """
    GPU-based whitening function that keeps the original input and output shapes.

    Args:
        X (torch.Tensor): Input tensor with shape (1, features, samples).
        rank (int): Number of components to keep.

    Returns:
        tuple: (whitening_matrix, whitened_X)
            - whitening_matrix is a tensor of shape (1, rank, features)
            - whitened_X is a tensor of shape (1, rank, samples)
    """
    # Center the data along the sample dimension.
    X_centered = X - X.mean(dim=2, keepdim=True)
    samples = X.shape[2]
    
    # Compute the covariance matrix of shape (1, features, features).
    cov = torch.matmul(X_centered, X_centered.transpose(1, 2)) / (samples - 1)
    
    # Compute the eigen decomposition of the covariance matrix.
    eigvals, eigvecs = torch.linalg.eigh(cov)
    
    # Rearrange eigenvalues and eigenvectors in descending order.
    eigvals = eigvals.flip(dims=[-1])
    eigvecs = eigvecs.flip(dims=[-1])
    
    # Stabilize eigenvector signs for reproducibility
    # Make the component with largest absolute value positive
    for i in range(eigvecs.shape[-1]):
        # Find index of maximum absolute value component for each batch and eigenvector
        max_idx = torch.argmax(torch.abs(eigvecs[:, :, i]), dim=1, keepdim=True)
        # Get the sign of that component
        signs = torch.sign(torch.gather(eigvecs[:, :, i], 1, max_idx))
        # Apply sign to make largest component positive
        eigvecs[:, :, i] *= signs
    
    # Select the top 'rank' components.
    eigvals = eigvals[:, :rank]
    eigvecs = eigvecs[:, :, :rank]
    
    # Compute scaling factors and form the whitening matrix.
    scale = 1.0 / torch.sqrt(eigvals + 1e-10)
    whitening_matrix = (eigvecs * scale.unsqueeze(1)).transpose(1, 2)
    
    # Compute the whitened data.
    whitened_X = torch.matmul(whitening_matrix, X_centered)
    
    return whitening_matrix, whitened_X


def find_ordering(S_list):
    """
    Find ordering of sources across views.
    Taken from https://github.com/hugorichard/multiviewica
    
    Args:
        S_list: List of source matrices
        
    Returns:
        Tuple of (u, orders, vals)
    """
    n_pb = len(S_list)
    
    # Normalize each source matrix
    for i in range(len(S_list)):
        S_list[i] /= np.linalg.norm(S_list[i], axis=1, keepdims=1)
    
    # Use first dataset as reference
    S = S_list[0].copy()
    n_components_ref = S.shape[0]
    
    # Store orderings for each dataset - list of arrays with different sizes
    orderings = [np.arange(n_components_ref)]  # First dataset keeps original order
    vals_list = []
    
    # Find best matching for other datasets
    for i, s in enumerate(S_list[1:]):
        # Compute correlation matrix between reference and current dataset
        # S shape: (n_components_ref, n_features)
        # s shape: (n_components_current, n_features)
        M = np.dot(S, s.T)  # Shape: (n_components_ref, n_components_current)
        
        # Find optimal assignment
        # This will match components from current dataset to reference
        ref_indices, current_indices = scipy.optimize.linear_sum_assignment(-abs(M))
        
        # Store the ordering for current dataset
        orderings.append(current_indices)
        
        # Store matching values
        vals = abs(M[ref_indices, current_indices])
        vals_list.append(vals)
    
    # Return the last assignment results for compatibility
    # But note that orderings now contains all the orderings
    return ref_indices, current_indices, vals if vals_list else (None, None, None)


def run_multi_view_ICA_on_datasets(
    datasets, a_values, c,
    batch_size=None,
    max_iter=10000,
    seed=0,
    mode='gpu',
    return_unmixing_matrices=False
):
    """
    Run multi-view ICA on multiple RNA-seq datasets.
    
    Args:
        datasets: List of DataFrames, each representing expression data for one view/species
                 Shape: (genes, samples)
        a_values: List of integers, number of components for each view/species
        c: Number of core components
        batch_size: Batch size for training (if None, uses full dataset)
        max_iter: Maximum iterations for optimizer
        seed: Random seed
        mode: 'gpu' or 'cpu'
        return_unmixing_matrices: If True, returns both sources and unmixing matrices
        
    Returns:
        If return_unmixing_matrices is False:
            List of DataFrames containing the source signals (samples x components) for each view
        If return_unmixing_matrices is True:
            Tuple of (sources, W_matrices, K_matrices) where:
            - sources: List of DataFrames (samples x components)
            - W_matrices: List of numpy arrays with unmixing matrices W
            - K_matrices: List of numpy arrays with whitening matrices K
    """
    # reproducibility - optimized seeding
    import os
    import random
    
    # Basic seeding for reproducibility
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        # Keep deterministic but allow benchmarking for better performance
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = True  # Allow cuDNN to find optimal algorithms
        # Enable TensorFloat32 for better performance on Ampere GPUs
        torch.set_float32_matmul_precision('high')

    # device
    device = torch.device("cuda:0" if mode == "gpu" and torch.cuda.is_available() else "cpu")
    
    # Set generator for device-specific operations
    if str(device) == "cuda:0":
        torch.cuda.set_device(device)
        # Clear GPU cache to prevent fragmentation
        torch.cuda.empty_cache()

    # to torch, transpose to (features, samples)
    X_np = [d.values.T.copy() for d in datasets]
    X_t = [torch.from_numpy(x).float().to(device) for x in X_np]

    # ranks
    ks = [torch.linalg.matrix_rank(x).item() for x in X_t]

    # add batch dim for whitening
    X_t = [x.unsqueeze(0) for x in X_t]
    
    # Parallel whitening if using GPU
    if str(device) == "cuda:0":
        # Process whitening in parallel on GPU
        with torch.cuda.stream(torch.cuda.Stream()):
            K, Xw = zip(*(whiten(x, k) for x, k in zip(X_t, ks)))
    else:
        K, Xw = zip(*(whiten(x, k) for x, k in zip(X_t, ks)))

    # (samples, rank)
    X_model = [xw[0].T.contiguous() for xw in Xw]

    # models with deterministic initialization
    models = []
    for i, (k, a) in enumerate(zip(ks, a_values)):
        # Set deterministic seed for this specific model
        model_seed = seed + i  # Simple deterministic seed per model
        torch.manual_seed(model_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(model_seed)
        
        model = MSIICA(k, a)
        
        # Explicitly reinitialize the orthogonal weights deterministically
        with torch.no_grad():
            # Create orthogonal matrix on CPU first for determinism
            orth_matrix = torch.nn.init.orthogonal_(torch.empty(a, k))
            model.W.weight.data = orth_matrix
        
        model = model.to(device)
        
        # Note: torch.compile removed due to recompilation issues with varying tensor sizes
        # The dynamic shapes cause cache misses and excessive recompilations
        
        models.append(model)

    # Use full dataset size if batch_size not specified
    if batch_size is None:
        batch_size = X_model[0].shape[0]  # number of samples

    # loader with optimizations
    generator = torch.Generator()
    generator.manual_seed(seed)
    
    # Note: X_model tensors are already on the device (GPU or CPU)
    # so we should NOT use pin_memory=True as it's only for CPU tensors
    train_loader = DataLoader(TensorDataset(*X_model),
                              batch_size=batch_size, 
                              shuffle=True, 
                              generator=generator,
                              pin_memory=False,  # Data is already on device
                              num_workers=0)  # Explicit for clarity

    # optimiser
    params = [p for m in models for p in m.parameters()]
    optimizer = torch.optim.LBFGS(params, line_search_fn="strong_wolfe",
                                  max_iter=max_iter)

    def train_closure(*inputs):
        def closure():
            optimizer.zero_grad()
            Sp = [m(x) for m, x in zip(models, inputs)]

            # ICA contrast
            loss = sum(torch.log(torch.cosh(s)).sum() for s in Sp)

            # anti-correlate shared part
            for i in range(len(Sp)):
                for j in range(i + 1, len(Sp)):
                    loss -= torch.trace(Sp[i][:, :c].T @ Sp[j][:, :c])

            loss.backward()

            # ensure contiguous grads for LBFGS
            for p in params:
                if p.grad is not None and not p.grad.is_contiguous():
                    p.grad = p.grad.contiguous()

            return loss
        return closure

    # training with GPU optimization
    for epoch in range(10):
        for batch in train_loader:
            # Data is already on the correct device from X_model
            optimizer.step(train_closure(*batch))
            
        # Periodic cache clearing to prevent memory fragmentation
        if str(device) == "cuda:0" and epoch % 3 == 0:
            torch.cuda.empty_cache()

    # final sources
    Sp_full = [m(x) for m, x in zip(models, X_model)]

    # (optional) reorder if you have a find_ordering()
    _ = find_ordering([s.detach().cpu().numpy().T for s in Sp_full])

    # to DataFrames
    dfs = [pd.DataFrame(s.detach().cpu().numpy()) for s in Sp_full]

    # variance-normalise columns
    for df in dfs:
        n = df.shape[0]
        centred = df - df.mean()
        var = centred.var(ddof=0).replace(0, 1)
        scale = np.sqrt(1.0 / (n * var))
        df.loc[:] = centred * scale
    
    if return_unmixing_matrices:
        # Extract W matrices from models and K matrices from whitening
        W_matrices = [m.W.weight.detach().cpu().numpy() for m in models]
        K_matrices = [k.detach().cpu().numpy() for k in K]
        return dfs, W_matrices, K_matrices
    else:
        return dfs


def run_multiview_ica(
    species_X_matrices: Dict[str, pd.DataFrame],
    a_values: Dict[str, int],
    c: int,
    mode: str = 'gpu',
    return_unmixing_matrices: bool = False,
    effect_size_threshold: Optional[float] = None,
    **kwargs
) -> Dict[str, pd.DataFrame]:
    """
    Run multi-view ICA natively using PyTorch.
    
    Args:
        species_X_matrices: Dictionary mapping species names to aligned X matrices
                           Each X matrix has shape (genes, samples)
        a_values: Dictionary mapping species names to number of components
        c: Number of core components
        mode: 'gpu' or 'cpu' mode
        return_unmixing_matrices: If True, returns unmixing matrices along with sources
        effect_size_threshold: Optional threshold for filtering components based on Cohen's d effect size.
                                 If provided, only keeps components above this threshold.
        **kwargs: Additional arguments passed to run_multi_view_ICA_on_datasets
        
    Returns:
        If return_unmixing_matrices is False:
            Dictionary mapping species names to M matrices (samples, components)
        If return_unmixing_matrices is True:
            Tuple of (M_matrices, W_matrices, K_matrices) where:
            - M_matrices: Dict mapping species to source matrices
            - W_matrices: Dict mapping species to unmixing matrices W
            - K_matrices: Dict mapping species to whitening matrices K
    """
    
    # Validate inputs
    n_species = len(species_X_matrices)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if len(a_values) != n_species:
        raise ValueError(f"Expected {n_species} a_values, got {len(a_values)}")
    
    # Prepare data with consistent ordering
    species_list = sorted(list(species_X_matrices.keys()))
    X_matrices = [species_X_matrices[sp] for sp in species_list]
    a_list = [a_values[sp] for sp in species_list]
    
    # Run multi-view ICA
    results = run_multi_view_ICA_on_datasets(
        X_matrices,
        a_list,
        c,
        mode=mode,
        return_unmixing_matrices=return_unmixing_matrices,
        **kwargs
    )
    
    if return_unmixing_matrices:
        # Unpack results
        result_dfs, W_list, K_list = results
        
        # Create dictionaries
        M_matrices = {}
        W_matrices = {}
        K_matrices = {}
        
        for species, df, W, K in zip(species_list, result_dfs, W_list, K_list):
            M_matrices[species] = df
            W_matrices[species] = W
            K_matrices[species] = K
        
        # Apply effect size filtering and formatting if threshold is provided
        if effect_size_threshold is not None:
            M_matrices, kept_core, kept_unique = _apply_effect_size_filtering_and_formatting(
                M_matrices, species_X_matrices, c, effect_size_threshold
            )
            print(f"Components saved: {kept_core} core, {kept_unique} unique")
        else:
            # Still apply formatting without filtering
            M_matrices = _apply_formatting_only(M_matrices, species_X_matrices, c)
        
        return M_matrices, W_matrices, K_matrices
    else:
        # Create results dictionary
        results_dict = {}
        for species, result in zip(species_list, results):
            results_dict[species] = result
        
        # Apply effect size filtering and formatting if threshold is provided
        if effect_size_threshold is not None:
            results_dict, kept_core, kept_unique = _apply_effect_size_filtering_and_formatting(
                results_dict, species_X_matrices, c, effect_size_threshold
            )
            print(f"Components saved: {kept_core} core, {kept_unique} unique")
        else:
            # Still apply formatting without filtering
            results_dict = _apply_formatting_only(results_dict, species_X_matrices, c)
        
        return results_dict


def _apply_effect_size_filtering_and_formatting(
    M_matrices: Dict[str, pd.DataFrame],
    species_X_matrices: Dict[str, pd.DataFrame], 
    c: int,
    effect_size_threshold: float
) -> Tuple[Dict[str, pd.DataFrame], int, int]:
    """
    Apply Cohen's d filtering and formatting to M matrices.
    
    Args:
        M_matrices: Dictionary mapping species to M matrices
        species_X_matrices: Dictionary mapping species to X matrices (for row indices)
        c: Number of core components
        effect_size_threshold: Cohen's d threshold for filtering components
        
    Returns:
        Tuple of (filtered_M_matrices, kept_core_count, kept_unique_count)
    """
    species_list = sorted(M_matrices.keys())
    
    # Calculate Cohen's d effect sizes for core components (average across species)
    core_effect_sizes = calculate_average_effect_sizes(
        {species: M_matrices[species].iloc[:, :c] for species in species_list}
    )
    
    # Filter core components
    core_keep_indices = [i for i, effect_size in enumerate(core_effect_sizes) if effect_size >= effect_size_threshold]
    kept_core = len(core_keep_indices)
    
    # Calculate Cohen's d effect sizes for unique components (individual per species)
    unique_keep_indices = {}
    total_kept_unique = 0
    
    for species in species_list:
        if M_matrices[species].shape[1] > c:  # Has unique components
            unique_components = M_matrices[species].iloc[:, c:]
            unique_effect_sizes = [
                calculate_cohens_d_effect_size(unique_components.iloc[:, i].values)
                for i in range(unique_components.shape[1])
            ]
            keep_indices = [i for i, effect_size in enumerate(unique_effect_sizes) if effect_size >= effect_size_threshold]
            unique_keep_indices[species] = keep_indices
            total_kept_unique += len(keep_indices)
        else:
            unique_keep_indices[species] = []
    
    # Filter and format M matrices
    filtered_M_matrices = {}
    for species in species_list:
        # Get row indices from corresponding X matrix
        row_indices = species_X_matrices[species].index
        
        # Filter core components
        if kept_core > 0:
            core_filtered = M_matrices[species].iloc[:, core_keep_indices]
            core_columns = [f"Core_{i+1}" for i in range(kept_core)]
        else:
            core_filtered = pd.DataFrame(index=row_indices)
            core_columns = []
        
        # Filter unique components for this species
        if unique_keep_indices[species]:
            unique_indices = [c + i for i in unique_keep_indices[species]]
            unique_filtered = M_matrices[species].iloc[:, unique_indices]
            unique_columns = [f"Unique_{i+1}" for i in range(len(unique_keep_indices[species]))]
        else:
            unique_filtered = pd.DataFrame(index=row_indices)
            unique_columns = []
        
        # Combine core and unique components
        if kept_core > 0 and unique_columns:
            filtered_df = pd.concat([core_filtered, unique_filtered], axis=1)
            column_names = core_columns + unique_columns
        elif kept_core > 0:
            filtered_df = core_filtered
            column_names = core_columns
        elif unique_columns:
            filtered_df = unique_filtered
            column_names = unique_columns
        else:
            filtered_df = pd.DataFrame(index=row_indices)
            column_names = []
        
        # Set proper row indices and column names
        filtered_df.index = row_indices
        filtered_df.columns = column_names
        
        filtered_M_matrices[species] = filtered_df
    
    return filtered_M_matrices, kept_core, total_kept_unique


def _apply_formatting_only(
    M_matrices: Dict[str, pd.DataFrame],
    species_X_matrices: Dict[str, pd.DataFrame],
    c: int
) -> Dict[str, pd.DataFrame]:
    """
    Apply formatting to M matrices without filtering.
    
    Args:
        M_matrices: Dictionary mapping species to M matrices
        species_X_matrices: Dictionary mapping species to X matrices (for row indices)
        c: Number of core components
        
    Returns:
        Dictionary of formatted M matrices
    """
    formatted_M_matrices = {}
    
    for species, M_matrix in M_matrices.items():
        # Get row indices from corresponding X matrix
        row_indices = species_X_matrices[species].index
        
        # Create column names
        total_components = M_matrix.shape[1]
        unique_components = total_components - c
        
        core_columns = [f"Core_{i+1}" for i in range(c)]
        unique_columns = [f"Unique_{i+1}" for i in range(unique_components)] if unique_components > 0 else []
        
        column_names = core_columns + unique_columns
        
        # Set proper row indices and column names
        formatted_df = M_matrix.copy()
        formatted_df.index = row_indices
        formatted_df.columns = column_names
        
        formatted_M_matrices[species] = formatted_df
    
    return formatted_M_matrices


