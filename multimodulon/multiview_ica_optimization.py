"""
Multi-view ICA optimization functions for selecting the optimal number of core components.

This module implements the NRE (Noise Reduction Error) method for automatically
determining the optimal number of core components in multi-view ICA.
"""

import os
import sys
import tempfile
import subprocess
import time
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Warning: matplotlib not available for plotting")


def create_modified_multiview_ica_script() -> str:
    """
    Create a modified version of the multiview ICA script that returns W and K matrices.
    
    Returns:
        str: The modified script content
    """
    script_content = '''
import torch
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
import torch.nn as nn
import geotorch
import scipy
import pickle

class MSIICA(nn.Module):
    def __init__(self, n_in, n_out, U=None, ortho=True):
        super().__init__()
        self.W = nn.Linear(n_in, n_out, bias=False)
        geotorch.orthogonal(self.W, "weight")
        if U is not None:
            self.W.weight = U.contiguous()
                
    def forward(self, Xw):
        S = self.W(Xw)
        return S

def find_ordering(S_list):
    # taken from https://github.com/hugorichard/multiviewica
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
            vals= abs(M[orders, u])

    return u,orders, vals

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

def run_multi_view_ICA_on_6_datasets_with_matrices(
    dataset_1_X, dataset_2_X, dataset_3_X, dataset_4_X, dataset_5_X, dataset_6_X,
    a1, a2, a3, a4, a5, a6, c,
    batch_size=None,
    max_iter=10000,
    seed=0,
    mode='gpu',
    return_whitening_unmixing=True
):
    """
    Modified version that also returns whitening and unmixing matrices.
    """
    # Set random seeds for reproducibility.
    torch.manual_seed(seed)
    np.random.seed(seed)

    # Set the device based on the mode input.
    if mode == 'gpu':
        device = torch.device("cuda:0")
    elif mode == 'cpu':
        device = torch.device("cpu")
    else:
        raise ValueError("Invalid mode. Choose 'gpu' or 'cpu'.")
    
    # Use full batch if not specified
    if batch_size is None:
        batch_size = len(dataset_1_X.index)
        
    # Convert DataFrames to numpy arrays and transpose to get shape (features, samples).
    X1_np = dataset_1_X.copy().values.T
    X2_np = dataset_2_X.copy().values.T
    X3_np = dataset_3_X.copy().values.T
    X4_np = dataset_4_X.copy().values.T
    X5_np = dataset_5_X.copy().values.T
    X6_np = dataset_6_X.copy().values.T
    
    # Store means for centering
    mu1 = np.mean(X1_np, axis=1, keepdims=True)
    mu2 = np.mean(X2_np, axis=1, keepdims=True)
    mu3 = np.mean(X3_np, axis=1, keepdims=True)
    mu4 = np.mean(X4_np, axis=1, keepdims=True)
    mu5 = np.mean(X5_np, axis=1, keepdims=True)
    mu6 = np.mean(X6_np, axis=1, keepdims=True)

    # Convert numpy arrays to torch tensors and move them to the device.
    X1_t = torch.from_numpy(X1_np).float().to(device)
    X2_t = torch.from_numpy(X2_np).float().to(device)
    X3_t = torch.from_numpy(X3_np).float().to(device)
    X4_t = torch.from_numpy(X4_np).float().to(device)
    X5_t = torch.from_numpy(X5_np).float().to(device)
    X6_t = torch.from_numpy(X6_np).float().to(device)

    # Compute the matrix rank for each domain using torch.
    k1 = torch.linalg.matrix_rank(X1_t).item()
    k2 = torch.linalg.matrix_rank(X2_t).item()
    k3 = torch.linalg.matrix_rank(X3_t).item()
    k4 = torch.linalg.matrix_rank(X4_t).item()
    k5 = torch.linalg.matrix_rank(X5_t).item()
    k6 = torch.linalg.matrix_rank(X6_t).item()

    # Reshape to (1, features, samples) for the whitening function.
    X1_t = X1_t.unsqueeze(0)
    X2_t = X2_t.unsqueeze(0)
    X3_t = X3_t.unsqueeze(0)
    X4_t = X4_t.unsqueeze(0)
    X5_t = X5_t.unsqueeze(0)
    X6_t = X6_t.unsqueeze(0)

    # Whiten each dataset.
    K1, Xw1 = whiten(X1_t, k1)
    K2, Xw2 = whiten(X2_t, k2)
    K3, Xw3 = whiten(X3_t, k3)
    K4, Xw4 = whiten(X4_t, k4)
    K5, Xw5 = whiten(X5_t, k5)
    K6, Xw6 = whiten(X6_t, k6)

    # Prepare model inputs: convert whitened data to shape (samples, rank).
    Xw1_model = Xw1[0].T.contiguous()
    Xw2_model = Xw2[0].T.contiguous()
    Xw3_model = Xw3[0].T.contiguous()
    Xw4_model = Xw4[0].T.contiguous()
    Xw5_model = Xw5[0].T.contiguous()
    Xw6_model = Xw6[0].T.contiguous()

    # Instantiate MSIICA models and move them to the device.
    model1 = MSIICA(k1, a1).to(device)
    model2 = MSIICA(k2, a2).to(device)
    model3 = MSIICA(k3, a3).to(device)
    model4 = MSIICA(k4, a4).to(device)
    model5 = MSIICA(k5, a5).to(device)
    model6 = MSIICA(k6, a6).to(device)

    # Create a DataLoader for the six datasets.
    train_data = TensorDataset(Xw1_model, Xw2_model, Xw3_model, Xw4_model, Xw5_model, Xw6_model)
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)

    # Combine model parameters for the LBFGS optimizer.
    params = (
        list(model1.parameters()) +
        list(model2.parameters()) +
        list(model3.parameters()) +
        list(model4.parameters()) +
        list(model5.parameters()) +
        list(model6.parameters())
    )
    optimizer = torch.optim.LBFGS(params, line_search_fn="strong_wolfe", max_iter=max_iter)

    def train_closure(x1, x2, x3, x4, x5, x6):
        def closure():
            optimizer.zero_grad()
            Sp1 = model1(x1)
            Sp2 = model2(x2)
            Sp3 = model3(x3)
            Sp4 = model4(x4)
            Sp5 = model5(x5)
            Sp6 = model6(x6)
    
            # Compute the loss
            loss = (
                torch.sum(torch.log(torch.cosh(Sp1))) +
                torch.sum(torch.log(torch.cosh(Sp2))) +
                torch.sum(torch.log(torch.cosh(Sp3))) +
                torch.sum(torch.log(torch.cosh(Sp4))) +
                torch.sum(torch.log(torch.cosh(Sp5))) +
                torch.sum(torch.log(torch.cosh(Sp6)))
            )
    
            # Subtract pairwise correlations
            Sp_list = [Sp1, Sp2, Sp3, Sp4, Sp5, Sp6]
            corr_val = 0.0
            for i in range(len(Sp_list)):
                for j in range(i + 1, len(Sp_list)):
                    corr_val += torch.trace(Sp_list[i][:, :c].T @ Sp_list[j][:, :c])
            loss = loss - corr_val
    
            loss.backward()
    
            # Ensure all gradients are contiguous
            for p in params:
                if p.grad is not None and not p.grad.is_contiguous():
                    p.grad.data = p.grad.data.contiguous()
    
            return loss
        return closure

    # Training loop.
    for _ in range(10):
        for batch in train_loader:
            optimizer.step(train_closure(*batch))

    # Obtain the final sources from each domain.
    Sp1_full = model1(Xw1_model)
    Sp2_full = model2(Xw2_model)
    Sp3_full = model3(Xw3_model)
    Sp4_full = model4(Xw4_model)
    Sp5_full = model5(Xw5_model)
    Sp6_full = model6(Xw6_model)

    # Optionally reorder sources
    S_list = [
        Sp1_full.detach().cpu().numpy().T,
        Sp2_full.detach().cpu().numpy().T,
        Sp3_full.detach().cpu().numpy().T,
        Sp4_full.detach().cpu().numpy().T,
        Sp5_full.detach().cpu().numpy().T,
        Sp6_full.detach().cpu().numpy().T
    ]
    _, _, _ = find_ordering(S_list)

    # Bring results back to the CPU and convert to pandas DataFrames.
    df1 = pd.DataFrame(Sp1_full.detach().cpu().numpy())
    df2 = pd.DataFrame(Sp2_full.detach().cpu().numpy())
    df3 = pd.DataFrame(Sp3_full.detach().cpu().numpy())
    df4 = pd.DataFrame(Sp4_full.detach().cpu().numpy())
    df5 = pd.DataFrame(Sp5_full.detach().cpu().numpy())
    df6 = pd.DataFrame(Sp6_full.detach().cpu().numpy())

    # Normalize each column to have variance 1/(number of rows)
    dataframes = [df1, df2, df3, df4, df5, df6]
    for df in dataframes:
        n_samples = df.shape[0]
        if n_samples == 0:
            continue
        # Center the data
        centered = df - df.mean()
        # Compute variance
        variances = centered.var(ddof=0)
        # Avoid division by zero (replace zero variances with 1)
        variances = variances.replace(0, 1)
        # Calculate scaling factor
        scaling_factors = np.sqrt(1 / (n_samples * variances))
        # Scale the data
        df_normalized = centered * scaling_factors
        # Update the DataFrame
        df.loc[:] = df_normalized.values

    if return_whitening_unmixing:
        # Extract unmixing matrices (W) from models
        W1 = model1.W.weight.detach().cpu().numpy()
        W2 = model2.W.weight.detach().cpu().numpy()
        W3 = model3.W.weight.detach().cpu().numpy()
        W4 = model4.W.weight.detach().cpu().numpy()
        W5 = model5.W.weight.detach().cpu().numpy()
        W6 = model6.W.weight.detach().cpu().numpy()
        
        # Convert whitening matrices to numpy
        K1_np = K1[0].detach().cpu().numpy()
        K2_np = K2[0].detach().cpu().numpy()
        K3_np = K3[0].detach().cpu().numpy()
        K4_np = K4[0].detach().cpu().numpy()
        K5_np = K5[0].detach().cpu().numpy()
        K6_np = K6[0].detach().cpu().numpy()
        
        # Store means
        means = [mu1, mu2, mu3, mu4, mu5, mu6]
        
        return (df1, df2, df3, df4, df5, df6, 
                [W1, W2, W3, W4, W5, W6],
                [K1_np, K2_np, K3_np, K4_np, K5_np, K6_np],
                means)
    else:
        return df1, df2, df3, df4, df5, df6
'''
    return script_content


def calculate_nre_proper(
    S_matrices: List[pd.DataFrame],
    k_core: int,
    normalize: bool = True
) -> float:
    """
    Calculate the proper Noise Reduction Error (NRE) for core components.
    
    The NRE measures how much the "core" components differ across views.
    Lower NRE means better core component sharing.
    
    Args:
        S_matrices: List of ICA source matrices (samples x components) for each view
        k_core: Number of core components to evaluate
        normalize: Whether to normalize by variance to avoid linear scaling with k
        
    Returns:
        NRE score (lower is better)
    """
    if k_core <= 0:
        return float('inf')
    
    n_views = len(S_matrices)
    n_samples = S_matrices[0].shape[0]
    
    # Extract first k_core components from each view
    core_components = []
    for S in S_matrices:
        if k_core > S.shape[1]:
            # If k_core exceeds available components, pad with zeros
            padded = np.zeros((n_samples, k_core))
            padded[:, :S.shape[1]] = S.values
            core_components.append(padded)
        else:
            core_components.append(S.values[:, :k_core])
    
    # Calculate NRE for each sample
    nre_per_sample = []
    
    for i in range(n_samples):
        # Extract core components for sample i across all views
        z_i_core_views = [comp[i, :] for comp in core_components]  # List of k_core arrays
        
        # Calculate mean across views for this sample
        z_bar_i = np.mean(z_i_core_views, axis=0)  # Shape: (k_core,)
        
        # Calculate residuals and sum of squared differences
        sum_squared_residuals = 0.0
        for z_i_d in z_i_core_views:
            residual = z_i_d - z_bar_i
            sum_squared_residuals += np.sum(residual ** 2)
        
        # Store raw residual sum for this sample
        nre_per_sample.append(sum_squared_residuals)
    
    # Calculate mean NRE across samples
    mean_nre = np.mean(nre_per_sample)
    
    if normalize:
        # Normalize by the total variance of core components to avoid linear scaling
        total_variance = 0.0
        for comp in core_components:
            # Calculate variance of each component and sum them
            comp_var = np.var(comp[:, :k_core], axis=0)
            total_variance += np.sum(comp_var)
        
        # Normalize by average variance per component per view
        avg_variance_per_comp = total_variance / (n_views * k_core) if k_core > 0 else 1.0
        
        # Normalized NRE: ratio of residual variance to total variance
        normalized_nre = mean_nre / (avg_variance_per_comp * k_core * n_views)
        return normalized_nre
    else:
        # Original NRE: divide by k_core
        return mean_nre / k_core


def calculate_alternative_metrics(
    S_matrices: List[pd.DataFrame],
    k_core: int
) -> Dict[str, float]:
    """
    Calculate alternative metrics for core component evaluation.
    
    Args:
        S_matrices: List of ICA source matrices for each view
        k_core: Number of core components to evaluate
        
    Returns:
        Dictionary with various metrics
    """
    if k_core <= 0:
        return {'correlation_score': 0.0, 'consistency_score': 0.0, 'stability_score': 0.0}
    
    n_views = len(S_matrices)
    
    # Extract core components
    core_components = []
    for S in S_matrices:
        if k_core > S.shape[1]:
            padded = np.zeros((S.shape[0], k_core))
            padded[:, :S.shape[1]] = S.values
            core_components.append(padded)
        else:
            core_components.append(S.values[:, :k_core])
    
    # 1. Correlation-based score: how well core components correlate across views
    correlation_scores = []
    for comp_idx in range(k_core):
        view_correlations = []
        for i in range(n_views):
            for j in range(i + 1, n_views):
                comp_i = core_components[i][:, comp_idx]
                comp_j = core_components[j][:, comp_idx]
                corr = np.abs(np.corrcoef(comp_i, comp_j)[0, 1])
                if not np.isnan(corr):
                    view_correlations.append(corr)
        
        if view_correlations:
            correlation_scores.append(np.mean(view_correlations))
    
    correlation_score = np.mean(correlation_scores) if correlation_scores else 0.0
    
    # 2. Consistency score: 1 - coefficient of variation across views
    consistency_scores = []
    for comp_idx in range(k_core):
        component_means = []
        for view_idx in range(n_views):
            comp_mean = np.mean(np.abs(core_components[view_idx][:, comp_idx]))
            component_means.append(comp_mean)
        
        if len(component_means) > 1 and np.mean(component_means) > 1e-8:
            cv = np.std(component_means) / np.mean(component_means)
            consistency_scores.append(max(0, 1 - cv))
        else:
            consistency_scores.append(0.0)
    
    consistency_score = np.mean(consistency_scores) if consistency_scores else 0.0
    
    # 3. Stability score: based on explained variance ratio
    total_variance = sum(np.sum(np.var(comp, axis=0)) for comp in core_components)
    core_variance = sum(np.sum(np.var(comp[:, :k_core], axis=0)) for comp in core_components)
    stability_score = core_variance / total_variance if total_variance > 0 else 0.0
    
    return {
        'correlation_score': correlation_score,
        'consistency_score': consistency_score,
        'stability_score': stability_score
    }


def run_nre_optimization_docker(
    species_X_matrices: Dict[str, pd.DataFrame],
    k_candidates: List[int],
    max_a_per_view: int,
    train_frac: float = 0.75,
    num_runs: int = 3,
    mode: str = 'gpu',
    docker_image: str = 'pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime',
    seed: int = 42
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Run NRE optimization using Docker to find optimal number of core components.
    
    Args:
        species_X_matrices: Dictionary mapping species names to X matrices
        k_candidates: List of candidate k values to try
        max_a_per_view: Maximum number of components per view
        train_frac: Fraction of data to use for training
        num_runs: Number of runs for averaging
        mode: 'gpu' or 'cpu'
        docker_image: Docker image to use
        seed: Random seed
        
    Returns:
        Tuple of (best_k, nre_dict, W_matrices_dict, K_matrices_dict, figure)
    """
    print("\n" + "="*60)
    print("NRE OPTIMIZATION FOR SHARED COMPONENTS")
    print("="*60)
    
    species_list = list(species_X_matrices.keys())
    if len(species_list) != 6:
        raise ValueError(f"Expected exactly 6 species, got {len(species_list)}")
    
    # Initialize storage for results
    mean_nre_per_k = {}
    all_nre_per_k = {k: [] for k in k_candidates}
    W_matrices_per_k = {}
    K_matrices_per_k = {}
    
    print(f"\nOptimization parameters:")
    print(f"  - k candidates: {k_candidates}")
    print(f"  - Train fraction: {train_frac}")
    print(f"  - Number of runs: {num_runs}")
    print(f"  - Max components per view: {max_a_per_view}")
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="nre_optimization_")
    
    # Save the modified script
    script_content = create_modified_multiview_ica_script()
    script_path = os.path.join(temp_dir, "multiview_ica_modified.py")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Create main execution script
    main_script = f'''
import sys
import os
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import train_test_split

sys.path.append('/workspace')
from multiview_ica_modified import run_multi_view_ICA_on_6_datasets_with_matrices

# Load data
species_list = {species_list}
X_views = []
for species in species_list:
    X = pd.read_csv(f"/workspace/data/{{species}}_X.csv", index_col=0)
    X_views.append(X)

k_candidates = {k_candidates}
train_frac = {train_frac}
num_runs = {num_runs}
max_a = {max_a_per_view}
mode = '{mode}'
seed = {seed}

results = {{}}

for run in range(num_runs):
    print(f"\\nRun {{run + 1}}/{{num_runs}}")
    
    # Split data with consistent indices across views
    n_samples = X_views[0].shape[0]
    indices = np.arange(n_samples)
    train_idx, test_idx = train_test_split(
        indices, train_size=train_frac, random_state=seed + run
    )
    
    X_train_views = [X.iloc[train_idx] for X in X_views]
    X_test_views = [X.iloc[test_idx] for X in X_views]
    
    for k in k_candidates:
        print(f"  Testing k={{k}}...")
        
        # Run multi-view ICA with current k
        results_ica = run_multi_view_ICA_on_6_datasets_with_matrices(
            *X_train_views,
            *[max_a] * 6,  # Use same a for all views
            k,  # c = k
            mode=mode,
            return_whitening_unmixing=True
        )
        
        # Extract results
        S_matrices = results_ica[:6]
        W_matrices = results_ica[6]
        K_matrices = results_ica[7]
        means_train = results_ica[8]
        
        # Calculate NRE on test data
        # Convert test data to numpy format expected by NRE calculation
        test_data_np = [df.values.T for df in X_test_views]
        
        # Calculate NRE manually in Python
        nre_scores_per_sample = []
        n_samples_test = X_test_views[0].shape[0]
        
        for i in range(n_samples_test):
            z_i_0_across_views = []
            
            for d in range(6):
                # Get the i-th sample
                x_id = X_test_views[d].iloc[i].values.reshape(-1, 1)
                
                # Center using training mean
                x_id_centered = x_id - means_train[d]
                
                # Apply whitening
                z_id_whitened = K_matrices[d] @ x_id_centered
                
                # Apply unmixing
                z_id = W_matrices[d] @ z_id_whitened.flatten()
                
                # Extract core components
                z_id_0 = z_id[:k]
                z_i_0_across_views.append(z_id_0)
            
            # Calculate mean across views
            z_i_0_stacked = np.stack(z_i_0_across_views)
            z_bar_i_0 = np.mean(z_i_0_stacked, axis=0)
            
            # Calculate squared differences
            sum_squared_diffs = 0.0
            for z_id_0 in z_i_0_across_views:
                z_hat_id = z_id_0 - z_bar_i_0
                sum_squared_diffs += np.sum(z_hat_id ** 2)
            
            nre_i = sum_squared_diffs / k
            nre_scores_per_sample.append(nre_i)
        
        nre = np.mean(nre_scores_per_sample)
        
        # Store results
        key = f"k_{{k}}_run_{{run}}"
        results[key] = {{
            'k': k,
            'run': run,
            'nre': nre,
            'W_matrices': W_matrices,
            'K_matrices': K_matrices,
            'means': means_train
        }}
        
        print(f"    NRE = {{nre:.6f}}")

# Save results
with open('/workspace/output/nre_results.pkl', 'wb') as f:
    pickle.dump(results, f)

print("\\nOptimization completed!")
'''
    
    main_script_path = os.path.join(temp_dir, "run_optimization.py")
    with open(main_script_path, 'w') as f:
        f.write(main_script)
    
    # Save X matrices
    data_dir = os.path.join(temp_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    for species, X in species_X_matrices.items():
        X.to_csv(os.path.join(data_dir, f"{species}_X.csv"))
    
    # Create output directory
    output_dir = os.path.join(temp_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Run Docker container
    print(f"\nRunning NRE optimization in Docker...")
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{temp_dir}:/workspace",
        "-w", "/workspace"
    ]
    
    if mode == 'gpu':
        docker_cmd.extend(["--gpus", "all"])
    
    docker_cmd.extend([
        docker_image,
        "python", "/workspace/run_optimization.py"
    ])
    
    start_time = time.time()
    result = subprocess.run(docker_cmd, capture_output=True, text=True, check=True)
    print(result.stdout)
    if result.stderr:
        print("Warnings:", result.stderr)
    
    # Load results
    with open(os.path.join(output_dir, "nre_results.pkl"), 'rb') as f:
        results = pickle.load(f)
    
    # Process results
    print("\nProcessing results...")
    for k in k_candidates:
        nre_values = []
        for run in range(num_runs):
            key = f"k_{k}_run_{run}"
            if key in results:
                nre_values.append(results[key]['nre'])
        
        if nre_values:
            mean_nre_per_k[k] = np.mean(nre_values)
            all_nre_per_k[k] = nre_values
            
            # Store W and K matrices from the last run
            last_key = f"k_{k}_run_{num_runs-1}"
            if last_key in results:
                W_matrices_per_k[k] = results[last_key]['W_matrices']
                K_matrices_per_k[k] = results[last_key]['K_matrices']
    
    # Find optimal k
    min_nre = min(mean_nre_per_k.values())
    # Find all k values within tolerance of minimum
    tolerance = 1e-8
    optimal_k_candidates = [k for k, nre in mean_nre_per_k.items() 
                           if abs(nre - min_nre) < tolerance]
    best_k = max(optimal_k_candidates)  # Choose largest k among ties
    
    print(f"\nOptimal k = {best_k} (NRE = {mean_nre_per_k[best_k]:.6f})")
    
    # Create plot
    if 'plt' in globals():
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        # Create a dummy figure object if matplotlib is not available
        class DummyFig:
            def savefig(self, *args, **kwargs): pass
        fig = DummyFig()
        ax = None
    
    # Plot mean NRE with error bars
    k_values = sorted(mean_nre_per_k.keys())
    mean_values = [mean_nre_per_k[k] for k in k_values]
    std_values = [np.std(all_nre_per_k[k]) if len(all_nre_per_k[k]) > 1 else 0 
                  for k in k_values]
    
    ax.errorbar(k_values, mean_values, yerr=std_values, 
                marker='o', markersize=8, capsize=5, capthick=2)
    
    # Highlight optimal k
    ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, 
               label=f'Optimal k = {best_k}')
    ax.scatter([best_k], [mean_nre_per_k[best_k]], color='red', s=100, zorder=5)
    
    ax.set_xlabel('Number of Shared Components (k)', fontsize=12)
    ax.set_ylabel('Mean NRE Score', fontsize=12)
    ax.set_title('NRE vs Number of Shared Components', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Clean up
    import shutil
    shutil.rmtree(temp_dir)
    
    total_time = time.time() - start_time
    print(f"\nTotal optimization time: {total_time:.2f}s")
    
    return best_k, mean_nre_per_k, W_matrices_per_k, K_matrices_per_k, fig


def run_nre_optimization_native(
    species_X_matrices: Dict[str, pd.DataFrame],
    k_candidates: List[int],
    max_a_per_view: int,
    train_frac: float = 0.75,
    num_runs: int = 3,
    mode: str = 'gpu',
    seed: int = 42
) -> Tuple[int, Dict[int, float], Dict[int, List], Dict[int, List], plt.Figure]:
    """
    Native implementation of NRE optimization using PyTorch.
    """
    # Import required functions - use the local implementation
    # Note: We avoid importing from the original multiview_ica due to syntax issues
    
    print("\n" + "="*60)
    print("NRE OPTIMIZATION FOR SHARED COMPONENTS (Native)")
    print("="*60)
    
    species_list = list(species_X_matrices.keys())
    n_species = len(species_list)
    if n_species < 2:
        raise ValueError(f"Expected at least 2 species, got {n_species}")
    
    print(f"Optimizing for {n_species} species: {species_list}")
    
    # Initialize storage for all metrics
    mean_nre_per_k = {}
    all_nre_per_k = {k: [] for k in k_candidates}
    mean_correlation_per_k = {}
    all_correlation_per_k = {k: [] for k in k_candidates}
    mean_consistency_per_k = {}
    all_consistency_per_k = {k: [] for k in k_candidates}
    mean_stability_per_k = {}
    all_stability_per_k = {k: [] for k in k_candidates}
    W_matrices_per_k = {}
    K_matrices_per_k = {}
    
    print(f"\nOptimization parameters:")
    print(f"  - k candidates: {k_candidates}")
    print(f"  - Train fraction: {train_frac}")
    print(f"  - Number of runs: {num_runs}")
    print(f"  - Max components per view: {max_a_per_view}")
    
    for run in range(num_runs):
        print(f"\nRun {run + 1}/{num_runs}")
        
        # Split data consistently across views
        n_samples = species_X_matrices[species_list[0]].shape[0]
        indices = np.arange(n_samples)
        train_idx, test_idx = train_test_split(
            indices, train_size=train_frac, random_state=seed + run
        )
        
        X_train_views = []
        X_test_views = []
        for species in species_list:
            X = species_X_matrices[species]
            X_train_views.append(X.iloc[train_idx])
            X_test_views.append(X.iloc[test_idx])
        
        for k in k_candidates:
            print(f"  Testing k={k}...")
            
            # Run multi-view ICA with current k
            # Import the generalized function
            from .multiview_ica import run_multi_view_ICA_on_datasets
            
            results = run_multi_view_ICA_on_datasets(
                X_train_views,
                [max_a_per_view] * n_species,  # Use same a for all views
                k,  # c = k (number of core components)
                mode=mode
            )
            
            # Train another model on test data to get comparable results
            test_results = run_multi_view_ICA_on_datasets(
                X_test_views,
                [max_a_per_view] * n_species,  # Use same a for all views
                k,  # c = k
                mode=mode
            )
            
            # Calculate proper NRE using the test results
            nre = calculate_nre_proper(test_results, k, normalize=True)
            
            # Also calculate alternative metrics for comparison
            alt_metrics = calculate_alternative_metrics(test_results, k)
            
            all_nre_per_k[k].append(nre)
            all_correlation_per_k[k].append(alt_metrics['correlation_score'])
            all_consistency_per_k[k].append(alt_metrics['consistency_score'])
            all_stability_per_k[k].append(alt_metrics['stability_score'])
            
            print(f"    NRE = {nre:.6f}")
            print(f"    Correlation = {alt_metrics['correlation_score']:.4f}")
            print(f"    Consistency = {alt_metrics['consistency_score']:.4f}")
            print(f"    Stability = {alt_metrics['stability_score']:.4f}")
    
    # Calculate mean metrics for each k
    for k in k_candidates:
        if all_nre_per_k[k]:
            mean_nre_per_k[k] = np.mean(all_nre_per_k[k])
            mean_correlation_per_k[k] = np.mean(all_correlation_per_k[k])
            mean_consistency_per_k[k] = np.mean(all_consistency_per_k[k])
            mean_stability_per_k[k] = np.mean(all_stability_per_k[k])
    
    # Find optimal k
    min_nre = min(mean_nre_per_k.values())
    tolerance = 1e-8
    optimal_k_candidates = [k for k, nre in mean_nre_per_k.items() 
                           if abs(nre - min_nre) < tolerance]
    best_k = max(optimal_k_candidates)
    
    print(f"\nOptimal k = {best_k} (NRE = {mean_nre_per_k[best_k]:.6f})")
    
    # Create comprehensive plot with all metrics
    if 'plt' in globals():
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
    else:
        # Create a dummy figure object if matplotlib is not available
        class DummyFig:
            def savefig(self, *args, **kwargs): pass
        fig = DummyFig()
        axes = [None, None, None, None]
    
    if axes[0] is not None:
        k_values = sorted(mean_nre_per_k.keys())
        
        # Prepare data for all metrics
        metrics_data = [
            {
                'name': 'NRE (Noise Reduction Error)',
                'mean_values': [mean_nre_per_k[k] for k in k_values],
                'std_values': [np.std(all_nre_per_k[k]) if len(all_nre_per_k[k]) > 1 else 0 for k in k_values],
                'color': 'blue',
                'ylabel': 'NRE Score',
                'optimal_value': mean_nre_per_k[best_k],
                'lower_is_better': True
            },
            {
                'name': 'Correlation Score',
                'mean_values': [mean_correlation_per_k[k] for k in k_values],
                'std_values': [np.std(all_correlation_per_k[k]) if len(all_correlation_per_k[k]) > 1 else 0 for k in k_values],
                'color': 'green',
                'ylabel': 'Correlation Score',
                'optimal_value': mean_correlation_per_k[best_k],
                'lower_is_better': False
            },
            {
                'name': 'Consistency Score',
                'mean_values': [mean_consistency_per_k[k] for k in k_values],
                'std_values': [np.std(all_consistency_per_k[k]) if len(all_consistency_per_k[k]) > 1 else 0 for k in k_values],
                'color': 'orange',
                'ylabel': 'Consistency Score',
                'optimal_value': mean_consistency_per_k[best_k],
                'lower_is_better': False
            },
            {
                'name': 'Stability Score',
                'mean_values': [mean_stability_per_k[k] for k in k_values],
                'std_values': [np.std(all_stability_per_k[k]) if len(all_stability_per_k[k]) > 1 else 0 for k in k_values],
                'color': 'purple',
                'ylabel': 'Stability Score',
                'optimal_value': mean_stability_per_k[best_k],
                'lower_is_better': False
            }
        ]
        
        # Plot each metric
        for i, (ax, metric) in enumerate(zip(axes, metrics_data)):
            # Plot the main curve with error bars
            ax.errorbar(k_values, metric['mean_values'], yerr=metric['std_values'], 
                       marker='o', markersize=6, capsize=4, capthick=1.5, 
                       color=metric['color'], alpha=0.8, linewidth=2)
            
            # Highlight optimal k
            ax.axvline(x=best_k, color='red', linestyle='--', alpha=0.7, linewidth=2)
            ax.scatter([best_k], [metric['optimal_value']], color='red', s=120, 
                      zorder=5, edgecolors='darkred', linewidth=2)
            
            # Formatting
            ax.set_xlabel('Number of Core Components (k)', fontsize=11)
            ax.set_ylabel(metric['ylabel'], fontsize=11)
            ax.set_title(f"{metric['name']}\n({'Lower is better' if metric['lower_is_better'] else 'Higher is better'})", 
                        fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # Add optimal value text
            direction = "↓" if metric['lower_is_better'] else "↑"
            ax.text(0.02, 0.98, f'Optimal k={best_k}: {metric["optimal_value"]:.4f} {direction}', 
                   transform=ax.transAxes, fontsize=10, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.8),
                   verticalalignment='top')
        
        # Add overall title
        fig.suptitle(f'Multi-view ICA Core Components Optimization\n'
                    f'Optimal k = {best_k} (NRE = {mean_nre_per_k[best_k]:.6f})', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        # Adjust layout
        plt.tight_layout(rect=[0, 0, 1, 0.94])
    
    return best_k, mean_nre_per_k, W_matrices_per_k, K_matrices_per_k, fig