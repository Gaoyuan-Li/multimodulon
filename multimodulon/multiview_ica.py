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
    Taken from https://github.com/hugorichard/multiviewica
    
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
        
    Returns:
        List of DataFrames containing the source signals (samples x components) for each view
    """
    # reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)

    # device
    device = torch.device("cuda:0" if mode == "gpu" and torch.cuda.is_available() else "cpu")

    # to torch, transpose to (features, samples)
    X_np = [d.values.T.copy() for d in datasets]
    X_t = [torch.from_numpy(x).float().to(device) for x in X_np]

    # ranks
    ks = [torch.linalg.matrix_rank(x).item() for x in X_t]

    # add batch dim for whitening
    X_t = [x.unsqueeze(0) for x in X_t]
    K, Xw = zip(*(whiten(x, k) for x, k in zip(X_t, ks)))

    # (samples, rank)
    X_model = [xw[0].T.contiguous() for xw in Xw]

    # models
    models = [MSIICA(k, a).to(device) for k, a in zip(ks, a_values)]

    # Use full dataset size if batch_size not specified
    if batch_size is None:
        batch_size = X_model[0].shape[0]  # number of samples

    # loader
    train_loader = DataLoader(TensorDataset(*X_model),
                              batch_size=batch_size, shuffle=True)

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

    # training
    for _ in range(10):
        for batch in train_loader:
            optimizer.step(train_closure(*batch))

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
        
    return dfs


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
                           Each X matrix has shape (genes, samples)
        a_values: Dictionary mapping species names to number of components
        c: Number of core components
        mode: 'gpu' or 'cpu' mode
        **kwargs: Additional arguments passed to run_multi_view_ICA_on_datasets
        
    Returns:
        Dictionary mapping species names to M matrices (samples, components)
    """
    
    # Validate inputs
    n_species = len(species_X_matrices)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    if len(a_values) != n_species:
        raise ValueError(f"Expected {n_species} a_values, got {len(a_values)}")
    
    # Prepare data
    species_list = list(species_X_matrices.keys())
    X_matrices = [species_X_matrices[sp] for sp in species_list]
    a_list = [a_values[sp] for sp in species_list]
    
    # Run multi-view ICA
    results = run_multi_view_ICA_on_datasets(
        X_matrices,
        a_list,
        c,
        mode=mode,
        **kwargs
    )
    
    # Create results dictionary
    results_dict = {}
    for species, result in zip(species_list, results):
        results_dict[species] = result
    
    return results_dict


# For backward compatibility, alias the docker function to native
run_multiview_ica_docker = run_multiview_ica_native