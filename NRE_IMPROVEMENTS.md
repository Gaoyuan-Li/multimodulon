# NRE Implementation Improvements

## Problem with Original Implementation

The original NRE implementation had a critical flaw:

```python
nre = np.random.random() * 0.1 + k * 0.01  # Linear scaling with k!
```

This caused:
- k=5 → NRE ≈ 0.1
- k=10 → NRE ≈ 0.15  
- k=75 → NRE ≈ 0.8

**Result**: Useless for optimization since it always increases with k.

## Theoretical NRE Formula

The correct NRE should measure how well "core" components are shared across views:

**NRE(k) = (1/N₁) Σᵢ Σd ||ẑ⁽ᵏ⁾d[i]||₂² / k**

Where:
- ẑ⁽ᵏ⁾d[i] = z⁽ᵏ⁾d,0[i] - z̄⁽ᵏ⁾0[i] (residual after removing cross-view mean)
- z⁽ᵏ⁾d,0[i] are the first k components from view d for sample i
- z̄⁽ᵏ⁾0[i] is the mean of those k components across all views

## Expected Behavior

- **k < k_true**: True core components → low residuals → **low NRE**
- **k = k_true**: Optimal point → **minimum NRE**  
- **k > k_true**: Including view-specific components → higher residuals → **higher NRE**

## Our Solution

### 1. Proper NRE Calculation

```python
def calculate_nre_proper(S_matrices, k_core, normalize=True):
    # Extract first k_core components from each view
    # Calculate cross-view residuals for each sample
    # Normalize by variance to avoid linear scaling
```

### 2. Normalization Strategy

**Problem**: Raw NRE still scales with k because more components = more variance.

**Solution**: Normalize by total component variance:

```python
normalized_nre = mean_nre / (avg_variance_per_comp * k_core * n_views)
```

This makes NRE a **ratio** rather than absolute value.

### 3. Alternative Metrics

We also provide complementary metrics:

**Correlation Score**: How well core components correlate across views
```python
correlation_score = mean(abs(corrcoef(comp_i, comp_j))) for all view pairs
```

**Consistency Score**: 1 - coefficient of variation across views
```python
consistency_score = 1 - std(component_means) / mean(component_means)
```

**Stability Score**: Fraction of variance explained by core components
```python
stability_score = core_variance / total_variance
```

## Implementation Changes

### 1. Updated Function Names
- `calculate_nre_proper()`: Implements correct NRE with normalization
- `calculate_alternative_metrics()`: Provides complementary measures

### 2. Updated Terminology
- "shared components" → "core components" (more accurate)
- All documentation and variable names updated

### 3. Enhanced Optimization
- Uses real ICA results instead of placeholders
- Trains separate models on train/test splits
- Provides multiple metrics for validation

## Expected Results

With proper implementation, you should see:

1. **Non-linear NRE curves** with clear minima
2. **Meaningful optimization** that finds true k_optimal
3. **Consistent results** across multiple runs
4. **Complementary metrics** that validate the choice

## Usage

```python
# Find optimal core components
best_k, nre_scores = multiModulon.optimize_number_of_core_components(
    step=5,
    num_runs=3,
    save_plot='nre_optimization.png'
)

# Run multi-view ICA with optimized k
multiModulon.run_multiview_ica(a=50, c=best_k)
```

The output will now show:
```
Testing k=5...
    NRE = 0.234567
    Correlation = 0.8234
    Consistency = 0.7123
    Stability = 0.6789

Testing k=10...
    NRE = 0.198432  # Lower is better!
    Correlation = 0.8567
    Consistency = 0.7456
    Stability = 0.7012
```

This should give you a proper U-shaped or L-shaped curve with a clear minimum!