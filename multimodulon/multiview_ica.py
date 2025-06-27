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
    
    Args:
        S_list: List of source matrices
        
    Returns:
        Tuple of (u, orders, vals)
    """
    n_pb = len(S_list)
    p = None
    for i in range(n_pb):
        p = S_list[i].shape[0] if p is None else np.min((p, S_list[i].shape[0]))

    for i in range(len(S_list)):
        S_list[i] /= np.linalg.norm(S_list[i], axis=1, keepdims=1)
    S = S_list[0].copy()
    order = np.arange(p)[None, :] * np.ones(n_pb, dtype=int)[:, None]
    for i, s in enumerate(S_list[1:]):
        M = np.dot(S, s.T)
        u, orders = scipy.optimize.linear_sum_assignment(-abs(M.T))
        order[i + 1] = orders
        vals = abs(M[orders, u])

    return u, orders, vals


def run_multi_view_ICA_on_datasets(
    datasets, a_values, c,
    batch_size=None,
    max_iter=10000,
    seed=0,
    mode='gpu'
):
    """
    Multi-view ICA implementation for any number of datasets.
    
    Args:
        datasets: List of DataFrames, each representing expression data for one view/species
        a_values: List of integers, number of components for each view/species
        c: Number of core components
        batch_size: Batch size for training
        max_iter: Maximum iterations for optimizer
        seed: Random seed
        mode: 'gpu' or 'cpu'
        
    Returns:
        List of DataFrames containing the ICA results for each view
    """
    # Validate inputs
    n_views = len(datasets)
    if len(a_values) != n_views:
        raise ValueError(f"Number of a_values ({len(a_values)}) must match number of datasets ({n_views})")
    
    print(f"Running multi-view ICA on {n_views} datasets")
    
    # Set random seeds for reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    # Set device
    if mode == 'gpu' and torch.cuda.is_available():
        device = torch.device("cuda:0")
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    else:
        device = torch.device("cpu")
        print("Using CPU")
    
    # Use full batch if not specified
    if batch_size is None:
        batch_size = len(datasets[0].index)
    
    # Convert DataFrames to numpy arrays and transpose to get shape (features, samples)
    X_tensors = []
    k_values = []
    K_matrices = []
    Xw_tensors = []
    Xw_models = []
    models = []
    
    for i, (dataset, a_val) in enumerate(zip(datasets, a_values)):
        print(f"Processing dataset {i+1}: {dataset.shape}")
        
        # Convert to numpy and transpose
        X_np = dataset.copy().values.T
        
        # Convert to torch tensor
        X_t = torch.from_numpy(X_np).float().to(device)
        X_tensors.append(X_t)
        
        # Compute matrix rank
        k = torch.linalg.matrix_rank(X_t).item()
        k_values.append(k)
        
        # Reshape for whitening
        X_t_reshaped = X_t.unsqueeze(0)
        
        # Whiten dataset
        K, Xw = whiten(X_t_reshaped, k)
        K_matrices.append(K)
        Xw_tensors.append(Xw)
        
        # Prepare model input
        Xw_model = Xw[0].T.contiguous()
        Xw_models.append(Xw_model)
        
        # Create model
        model = MSIICA(k, a_val).to(device)
        models.append(model)
    
    # Create DataLoader
    train_data = TensorDataset(*Xw_models)
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    
    # Set up optimizer - collect all model parameters
    params = []
    for model in models:
        params.extend(list(model.parameters()))
    optimizer = torch.optim.LBFGS(params, line_search_fn="strong_wolfe", max_iter=max_iter)
    
    def train_closure(*batch_data):
        def closure():
            optimizer.zero_grad()
            
            # Forward pass through all models
            Sp_list = []
            for i, (model, x_batch) in enumerate(zip(models, batch_data)):
                Sp = model(x_batch)
                Sp_list.append(Sp)
            
            # Compute loss - sum of log cosh for all views
            loss = 0.0
            for Sp in Sp_list:
                loss += torch.sum(torch.log(torch.cosh(Sp)))
            
            # Subtract pairwise correlations for core components
            corr_val = 0.0
            for i in range(len(Sp_list)):
                for j in range(i + 1, len(Sp_list)):
                    # Use only the first c components for correlation
                    corr_val += torch.trace(Sp_list[i][:, :c].T @ Sp_list[j][:, :c])
            loss = loss - corr_val
            
            loss.backward()
            
            # Ensure gradients are contiguous
            for p in params:
                if p.grad is not None and not p.grad.is_contiguous():
                    p.grad.data = p.grad.data.contiguous()
            
            return loss
        return closure
    
    # Training loop
    print("Starting training...")
    for epoch in range(10):
        for batch in train_loader:
            optimizer.step(train_closure(*batch))
        if epoch % 2 == 0:
            print(f"Epoch {epoch}/10 completed")
    
    # Get final results
    print("Computing final results...")
    results = []
    for i, (model, Xw_model) in enumerate(zip(models, Xw_models)):
        Sp_full = model(Xw_model)
        df = pd.DataFrame(Sp_full.detach().cpu().numpy())
        results.append(df)
        print(f"Dataset {i+1} result shape: {df.shape}")
    
    # Normalize DataFrames
    for df in results:
        n_samples = df.shape[0]
        if n_samples == 0:
            continue
        centered = df - df.mean()
        variances = centered.var(ddof=0)
        variances = variances.replace(0, 1)
        scaling_factors = np.sqrt(1 / (n_samples * variances))
        df_normalized = centered * scaling_factors
        df.loc[:] = df_normalized.values
    
    return results


def run_multiview_ica_native(
    species_X_matrices: Dict[str, pd.DataFrame],
    a_values: Dict[str, int],
    c: int,
    mode: str = 'gpu',
    **kwargs
) -> Dict[str, pd.DataFrame]:
    """
    Run multi-view ICA natively using PyTorch.
    
    Args:
        species_X_matrices: Dictionary mapping species names to aligned X matrices
        a_values: Dictionary mapping species names to number of components
        c: Number of core components
        mode: 'gpu' or 'cpu' mode
        
    Returns:
        Dictionary mapping species names to M matrices (ICA results)
    """
    print("\n" + "="*60)
    print("MULTI-VIEW ICA EXECUTION (Native)")
    print("="*60)
    
    start_time = time.time()
    
    # Validate inputs
    n_species = len(species_X_matrices)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if len(a_values) != n_species:
        raise ValueError(f"Expected {n_species} a_values, got {len(a_values)}")
    
    # Prepare data
    print("\n[1/2] Preparing data...")
    prep_start = time.time()
    
    species_list = list(species_X_matrices.keys())
    X_matrices = [species_X_matrices[sp] for sp in species_list]
    a_list = [a_values[sp] for sp in species_list]
    
    print(f"Species: {species_list}")
    print(f"a_values: {a_list}")
    print(f"c (core components): {c}")
    print(f"Mode: {mode}")
    print(f"✓ Data prepared in {time.time() - prep_start:.2f}s")
    
    # Run multi-view ICA
    print(f"\n[2/2] Running multi-view ICA...")
    run_start = time.time()
    
    results = run_multi_view_ICA_on_datasets(
        X_matrices,
        a_list,
        c,
        mode=mode,
        **kwargs
    )
    
    print(f"✓ Multi-view ICA completed in {time.time() - run_start:.2f}s")
    
    # Create results dictionary
    results_dict = {}
    for species, result in zip(species_list, results):
        results_dict[species] = result
        print(f"✓ {species} M matrix: {result.shape}")
    
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"TOTAL EXECUTION TIME: {total_time:.2f}s")
    print(f"{'='*60}\n")
    
    return results_dict


# For backward compatibility, alias the docker function to native
run_multiview_ica_docker = run_multiview_ica_native