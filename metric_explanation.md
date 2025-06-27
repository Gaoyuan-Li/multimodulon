# Multi-view ICA Optimization: NRE Metric Explanation

This document provides detailed explanation of the NRE (Normalized Reconstruction Error) metric used for optimizing the number of core components in multi-view ICA analysis.

## Overview

When running `optimize_number_of_core_components()`, the algorithm uses the **NRE (Normalized Reconstruction Error)** metric as defined in the multi-view ICA paper to determine the optimal number of core components.

## NRE Definition

### NRE (Normalized Reconstruction Error) 🔴 Lower is Better

**Purpose**: Measures how well the first k components are truly "core" (shared) across all views.

**Mathematical Definition**: As defined in the multi-view ICA paper:

**NRE(k) := (1/N₁) Σᵢ₌₁^{N₁} Σ_{d=1}^D ||ẑ_d^{(k)}_i||² / k**

Where:
- **ẑ_d^{(k)}_i = z_{d,0}^{(k)}_i - z̄_0^{(k)}_i** (residual for sample i after removing cross-view mean)
- **z_{d,0}^{(k)}_i** are the first k components from view d for sample i
- **z̄_0^{(k)}_i = (1/D) Σ_{ℓ=1}^D z_{ℓ,0}^{(k)}_i** is the mean across views for sample i
- **D** is the number of views (species)
- **N₁** is the number of test samples

**Implementation**:
```python
def calculate_nre_proper(S_matrices, k_core):
    D = len(S_matrices)  # Number of views
    N_1 = S_matrices[0].shape[0]  # Number of test samples
    
    # Extract first k_core components from each view (z_{d,0}^(k))
    z_d_0_k = []
    for S in S_matrices:
        z_d_0_k.append(S.values[:, :k_core])
    
    # Calculate NRE for each sample individually, then average
    nre_per_sample = []
    
    for i in range(N_1):
        # Extract components for sample i from all views
        z_i_views = [z_d[i, :] for z_d in z_d_0_k]
        
        # Calculate bar_z_0^(k)_i = (1/D) * sum_{l=1}^D z_{l,0}^(k)_i
        bar_z_i = np.mean(z_i_views, axis=0)
        
        # Calculate NRE for this sample
        sample_nre = 0.0
        for z_i_d in z_i_views:
            hat_z_i_d = z_i_d - bar_z_i
            sample_nre += np.sum(hat_z_i_d**2)
        
        nre_per_sample.append(sample_nre / k_core)
    
    # Average over all test samples
    return np.mean(nre_per_sample)
```

**Expected Behavior**:
- k < k_optimal: Low NRE (true core components are well-shared)
- k = k_optimal: Minimum NRE (optimal balance)
- k > k_optimal: Higher NRE (including view-specific components increases residuals)

**Interpretation**: Look for the k value where NRE reaches its minimum in the optimization plot.

## Theoretical Background

According to the multi-view ICA paper, the NRE metric is designed to evaluate the quality of shared source estimates for varying values of k. The key insight is:

- If **k = k*** (the true number of core components), the estimated core components **z_{d,0}^{(k*)}** become well-aligned across views
- This means **ẑ_d^{(k*)}** (the residuals) follow a normal distribution with mean 0 and variance **(D-1)/D · σ² · I_{k*}**
- The NRE score approaches a theoretical minimum at the optimal k

### Selection Strategy

As specified in the paper, we use the **elbow method** with the following criterion:

**k* = max{argmin_k NRE(k)}**

This means we choose the **largest k** among all values that achieve the minimum NRE score, accounting for numerical tolerance.

## Cross-Validation Procedure

The optimization follows the paper's methodology:

1. **Data Split**: Divide each view's data into training (75%) and test (25%) sets
2. **Model Training**: Train multi-view ICA models on training data for each candidate k
3. **NRE Evaluation**: Calculate NRE using the test data to avoid overfitting
4. **Multiple Runs**: Repeat the process multiple times with different random splits for robustness

## How to Interpret the Optimization Plot

The optimization plot shows:
- **X-axis**: Number of core components (k)
- **Y-axis**: NRE Score (lower is better)
- **Error bars**: Standard deviation across multiple cross-validation runs
- **Red vertical line**: Optimal k value
- **Red dot**: NRE value at optimal k

### Typical Patterns:

**Good Optimization**:
- NRE shows a clear **U-shape** or **L-shape** with a distinct minimum
- The minimum occurs at a reasonable k value (not at the extremes)
- Error bars are relatively small, indicating consistent results

**Poor Optimization** (may need parameter adjustment):
- NRE continuously decreases without a clear minimum
- Very high error bars indicating unstable results
- Minimum occurs at k=1 or at the maximum tested k

## Example Interpretation

```
Testing k=5...
    NRE = 0.234567

Testing k=10...
    NRE = 0.198432        # Lower - potential optimum

Testing k=15...
    NRE = 0.245123        # Higher - past optimum

Optimal k = 10 (NRE = 0.198432)
```

In this example, **k=10** is optimal because NRE reaches its minimum there.

## Advanced Usage Tips

1. **Step Size**: Use smaller steps (e.g., `step=2`) for fine-tuning around the optimal region
2. **Number of Runs**: Increase `num_runs` for more reliable optimization with noisy data
3. **Training Fraction**: Adjust `train_frac` based on data size (more training data for small datasets)
4. **Visual Inspection**: Always examine the plot - the automated selection may miss nuances

## Implementation Usage

```python
# Find optimal core components
best_k, nre_scores = multiModulon.optimize_number_of_core_components(
    step=5,              # Test k = 5, 10, 15, 20, ...
    num_runs=3,          # 3 cross-validation runs
    train_frac=0.75,     # 75% train, 25% test
    save_plot='nre_optimization.png'
)

print(f"Optimal number of core components: {best_k}")

# Run multi-view ICA with optimized k
multiModulon.run_multiview_ica(a=50, c=best_k)
```

## References

- Multi-view ICA paper methodology for NRE calculation
- Cross-validation best practices for model selection
- Elbow method for optimization criterion selection