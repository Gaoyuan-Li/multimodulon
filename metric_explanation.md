# Multi-view ICA Optimization Metrics Explanation

This document provides detailed explanations of the four metrics used for optimizing the number of core components in multi-view ICA analysis.

## Overview

When running `optimize_number_of_core_components()`, the algorithm evaluates four complementary metrics to determine the optimal number of core components:

1. **NRE (Noise Reduction Error)** - Primary optimization metric
2. **Correlation Score** - Cross-view component similarity 
3. **Consistency Score** - Stability of component values across views
4. **Stability Score** - Fraction of variance captured by core components

## Metric Definitions

### 1. NRE (Noise Reduction Error) 🔴 Lower is Better

**Purpose**: Measures how well the first k components are truly "core" (shared) across all views.

**Theory**: True core components should have similar patterns across different species/views. When we include too many components (k > optimal), we start including view-specific components that don't share patterns, increasing the residual error.

**Implementation**:
```python
def calculate_nre_proper(S_matrices, k_core, normalize=True):
    residuals = []
    for sample_idx in range(n_samples):
        # Extract first k components for this sample across all views
        sample_components = [S[sample_idx, :k_core] for S in S_matrices]
        
        # Calculate cross-view mean for this sample
        cross_view_mean = np.mean(sample_components, axis=0)
        
        # Calculate residuals (how much each view deviates from the mean)
        for comp in sample_components:
            residual = comp - cross_view_mean
            residuals.append(np.sum(residual**2))
    
    mean_nre = np.mean(residuals)
    
    if normalize:
        # Normalize by component variance to avoid linear scaling with k
        avg_variance = np.mean([np.var(S[:, :k_core]) for S in S_matrices])
        normalized_nre = mean_nre / (avg_variance * k_core * len(S_matrices))
        return normalized_nre
    
    return mean_nre
```

**Expected Behavior**:
- k < k_optimal: Low NRE (true core components are well-shared)
- k = k_optimal: Minimum NRE (optimal balance)
- k > k_optimal: Higher NRE (including view-specific components increases residuals)

**Interpretation**: Look for the k value where NRE reaches its minimum in the optimization plot.

---

### 2. Correlation Score 🟢 Higher is Better

**Purpose**: Measures how well corresponding core components correlate between different views.

**Theory**: True core components should show strong correlations when comparing the same component index across different species/views.

**Implementation**:
```python
def calculate_correlation_score(S_matrices, k_core):
    correlations = []
    
    # For each component index in the core set
    for comp_idx in range(k_core):
        # For each pair of views
        for i in range(len(S_matrices)):
            for j in range(i + 1, len(S_matrices)):
                # Extract the same component from both views
                comp_i = S_matrices[i][:, comp_idx]
                comp_j = S_matrices[j][:, comp_idx]
                
                # Calculate correlation coefficient
                corr = np.corrcoef(comp_i, comp_j)[0, 1]
                correlations.append(abs(corr))  # Use absolute value
    
    return np.mean(correlations)
```

**Expected Behavior**:
- k < k_optimal: High correlation (all components are truly core)
- k = k_optimal: High correlation maintained
- k > k_optimal: Decreasing correlation (including view-specific components)

**Interpretation**: Look for k values that maintain high correlation scores (> 0.7-0.8).

---

### 3. Consistency Score 🟢 Higher is Better

**Purpose**: Measures how consistent component magnitudes are across different views.

**Theory**: Core components should have similar magnitudes/scales across views. Large variations in magnitude suggest view-specific scaling or non-core components.

**Implementation**:
```python
def calculate_consistency_score(S_matrices, k_core):
    consistency_scores = []
    
    # For each component in the core set
    for comp_idx in range(k_core):
        # Calculate mean magnitude for this component in each view
        component_means = []
        for S in S_matrices:
            comp_mean = np.mean(np.abs(S[:, comp_idx]))
            component_means.append(comp_mean)
        
        # Calculate coefficient of variation (CV = std/mean)
        mean_of_means = np.mean(component_means)
        std_of_means = np.std(component_means)
        
        if mean_of_means > 0:
            cv = std_of_means / mean_of_means
            consistency = 1 - cv  # Higher consistency = lower CV
        else:
            consistency = 0
        
        consistency_scores.append(max(0, consistency))  # Clamp to [0,1]
    
    return np.mean(consistency_scores)
```

**Expected Behavior**:
- k < k_optimal: High consistency (core components have similar scales)
- k = k_optimal: High consistency maintained
- k > k_optimal: Decreasing consistency (view-specific components vary in scale)

**Interpretation**: Look for k values where consistency remains high (> 0.6-0.8).

---

### 4. Stability Score 🟢 Higher is Better

**Purpose**: Measures what fraction of the total variance is captured by the core components.

**Theory**: Core components should explain a substantial portion of the variance. If the score is too low, we might be missing important shared patterns.

**Implementation**:
```python
def calculate_stability_score(S_matrices, k_core):
    stability_scores = []
    
    for S in S_matrices:
        # Calculate variance of core components (first k)
        core_variance = np.sum(np.var(S[:, :k_core], axis=0))
        
        # Calculate total variance across all components
        total_variance = np.sum(np.var(S, axis=0))
        
        if total_variance > 0:
            stability = core_variance / total_variance
        else:
            stability = 0
        
        stability_scores.append(stability)
    
    return np.mean(stability_scores)
```

**Expected Behavior**:
- k increases: Stability score increases (more variance captured)
- k = k_optimal: Good balance between capturing variance and avoiding overfitting
- k > k_optimal: Diminishing returns (additional view-specific components don't add much shared variance)

**Interpretation**: Look for k values where stability reaches a reasonable plateau (> 0.5-0.7).

---

## How to Interpret the Optimization Plot

The optimization plot shows all four metrics in a 2×2 grid:

```
NRE (Lower is Better)     |  Correlation Score (Higher is Better)
--------------------------|----------------------------------
Consistency Score         |  Stability Score (Higher is Better)
(Higher is Better)        |
```

### Finding the Optimal k:

1. **Primary**: Look for the k value where **NRE reaches its minimum** (red line and marker)
2. **Validation**: Check that at this k value:
   - **Correlation** remains high (> 0.7)
   - **Consistency** remains high (> 0.6) 
   - **Stability** shows reasonable variance capture (> 0.5)

### Typical Patterns:

**Good Optimization**:
- NRE shows a clear U-shape or L-shape with a minimum
- Correlation stays high until k becomes too large
- Consistency remains stable until k becomes too large
- Stability increases smoothly and plateaus

**Poor Optimization** (may need different parameters):
- NRE continuously decreases (no clear minimum)
- All metrics show monotonic trends without clear optima
- Very low correlation scores (< 0.5) across all k values

## Example Interpretation

```
Testing k=5...
    NRE = 0.234567        # Moderate
    Correlation = 0.8234  # Good - components correlate well across views
    Consistency = 0.7123  # Good - similar magnitudes across views  
    Stability = 0.4789    # Moderate - capturing ~48% of variance

Testing k=10...
    NRE = 0.198432        # Lower is better - this might be optimal!
    Correlation = 0.8567  # Still high - good sign
    Consistency = 0.7456  # Still high - good sign
    Stability = 0.6012    # Higher - capturing more variance

Testing k=15...
    NRE = 0.245123        # Increased - probably too many components
    Correlation = 0.7234  # Decreased - including view-specific components
    Consistency = 0.6789  # Decreased - more variation in scales
    Stability = 0.7234    # Higher but diminishing returns
```

In this example, **k=10** would be optimal: NRE is minimized while maintaining high correlation and consistency.

## Advanced Usage Tips

1. **Cross-validation**: The algorithm runs multiple train/test splits to ensure robust results
2. **Step size**: Use smaller steps (e.g., step=2) for fine-tuning around the optimal region
3. **Multiple runs**: Increase `num_runs` for more reliable optimization in noisy data
4. **Manual validation**: Always inspect the plot visually - automated selection may miss nuances

## References

- NRE methodology adapted from multi-view ICA literature
- Correlation and consistency metrics based on component analysis best practices
- Stability scoring follows standard variance decomposition principles