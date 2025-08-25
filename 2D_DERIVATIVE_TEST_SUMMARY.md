# 2D Derivative Matrix Testing - Results Summary

## Overview
This document summarizes the comprehensive testing of 2D derivative matrices constructed using tensor products of 1D Legendre-Gauss-Lobatto (LGL) derivative matrices.

## Test Setup

### Implementation Method
- **1D Base**: Uses `derv()` subroutine from `lgl_baseline.f90` to create 1D differentiation matrix
- **2D Extension**: Tensor product construction:
  - `D_x = D_1D ⊗ I_y` (differentiate in x, identity in y)
  - `D_y = I_x ⊗ D_1D` (identity in x, differentiate in y)
- **Grid Mapping**: `node_index = i + j*(N+1)` for 2D tensor product grid

### Test Functions
1. **Trigonometric**: `f(x,y) = sin(πx)cos(πy)`
2. **Polynomial**: `f(x,y) = x³y² + xy³` (lower order), `f(x,y) = x⁴y³ + x²y⁵` (higher order)
3. **Exponential**: `f(x,y) = exp(x+y)` and `f(x,y) = exp(0.5*(x+y))`

## Key Results

### Spectral Convergence Achieved
Both even and odd polynomial orders demonstrate **exponential convergence** for smooth functions:

#### Even Orders (N = 2, 4, 6, ..., 20)
```
Order N   Max Error dx   Max Error dy   DOF
-------   ------------   ------------   ---
   2      3.14E+00       4.90E-16        9
   4      1.59E+00       8.67E-01       25
   8      1.08E-02       3.25E-03       81
  12      7.46E-06       1.72E-06      169
  20      1.18E-12       3.45E-12      441
```

#### Odd Orders (N = 3, 5, 7, ..., 21)
```
Order N   Max Error dx   Max Error dy   DOF
-------   ------------   ------------   ---
   3      2.37E+00       2.87E+00       16
   5      3.14E-01       4.86E-01       36
   9      6.15E-04       2.06E-03      100
  13      2.04E-07       9.40E-07      196
  21      3.02E-14       1.10E-13      484
```

### Convergence Rates
- **Trigonometric functions**: Rates of 10-40 per doubling of polynomial order
- **Polynomial functions**: Machine precision achieved when N ≥ polynomial degree
- **Exponential functions**: Consistent exponential convergence

## Key Differences: Even vs Odd Orders

### Mathematical Properties

#### Even Orders (N = 2k)
- **Grid points**: Always include `x = 0, y = 0`
- **Symmetry**: Natural for symmetric functions
- **Example N=4**: `x_1d = [-1.0, -0.6547, 0.0, +0.6547, +1.0]`

#### Odd Orders (N = 2k+1)  
- **Grid points**: NO point at `x = 0, y = 0`
- **Symmetry**: Points symmetric about origin but no center point
- **Example N=5**: `x_1d = [-1.0, -0.7651, -0.2852, +0.2852, +0.7651, +1.0]`

### Accuracy Comparison
For similar computational cost (adjacent orders):
- **Even orders**: Often slightly more accurate for symmetric functions
- **Odd orders**: More efficient (better accuracy per DOF) for many cases
- **DOF Ratio**: Odd/Even ≈ 1.1-1.4 for adjacent orders

### Performance Metrics
```
Comparison   DOF Ratio   Error Ratio dx   Error Ratio dy
4/5 points      1.44         5.05            1.78
8/9 points      1.23        17.53            1.58
12/13 points    1.16        36.55            1.83
```

## LGL Collocation Points

### Mathematical Definition
Points are roots of `(1-x²)P'_N(x) = 0` where `P_N(x)` is the N-th Legendre polynomial.

### Key Properties
- **Boundary conditions**: Always include `x = ±1`
- **Clustering**: Points cluster near boundaries for better resolution
- **High precision**: Newton-Raphson with `ε = 10⁻¹²` convergence
- **Symmetry**: `x_i = -x_{N-i}` (mirror symmetry)

### Grid Point Distribution Examples
```
N=2:  [-1.000,  0.000,  1.000]
N=5:  [-1.000, -0.765, -0.285,  0.285,  0.765,  1.000]
N=8:  [-1.000, -0.900, -0.677, -0.363,  0.000,  0.363,  0.677,  0.900,  1.000]
```

## Implementation Validation

### Test Results Summary
✅ **Polynomial exactness**: Machine precision for polynomials of degree ≤ N  
✅ **Smooth function accuracy**: Exponential convergence demonstrated  
✅ **Tensor product structure**: Correct 2D derivative matrix construction  
✅ **Spectral accuracy**: Expected high-order convergence rates achieved  

### Memory and Computational Complexity
- **Matrix size**: `(N+1)² × (N+1)²` for each derivative direction
- **Sparsity**: Each row has exactly `(N+1)` non-zero entries
- **Storage**: Efficient for N ≤ 20, consider sparse storage for higher orders

## Recommendations

### When to Use Even vs Odd Orders
- **Even orders**: Preferred for symmetric problems, functions with behavior at origin
- **Odd orders**: Often more efficient, good for general applications
- **High accuracy**: Both achieve spectral accuracy; choice depends on specific problem

### Practical Guidelines
- **N = 8-12**: Good balance of accuracy vs computational cost
- **N > 16**: Consider numerical stability and computational resources
- **Memory**: 2D matrices become large; consider matrix-free implementations for very high orders

## Code Files Created
1. `test_derv_2D.f90` - Basic 2D derivative matrix testing
2. `test_derv_2D_enhanced.f90` - Multiple test functions and convergence analysis
3. `test_derv_2D_odd_orders.f90` - Comprehensive odd order testing
4. `compare_even_odd_orders.f90` - Direct comparison of even vs odd orders
5. `show_x1d_points.f90` - Display LGL collocation points

All tests confirm that the tensor product approach using the 1D `derv()` subroutine from `lgl_baseline.f90` successfully creates accurate 2D derivative matrices with spectral convergence properties.
