# Data dictionary

The values describe normalized mathematical models, not physical measurements. JSON files record floating-point numbers and parameters. Plot files are whitespace-separated text with header rows. The source programs provide the exact computational formulas and assertion tolerances.

## Interval results

`verification/continuum_results.json` records the domain, coupling, trial basis, software versions and two collections. `cases` contains three finest-resolution parameter cases with twelve collar-width rows each. `refinement` contains twelve fractional-order/degree/quadrature configurations, including the three finest cases. `degree` is the largest Legendre degree; `quadrature_order` is the outer quadrature order and `interior_x_quadrature_order` the inner moment order. `S_galerkin` is the computed finite-space corrector energy; `G_exact` is 1/s.

| Row field | Meaning |
|---|---|
| `epsilon` | Exterior collar width. |
| `rho` | epsilon^(1-2s). |
| `lambda1_galerkin` | Principal eigenvalue of the quadrature-assembled Galerkin matrix, evaluated by the scalar Schur equation. |
| `constant_trial` | Collar energy at the normalized constant trial function. |
| `trial_minus_lambda` | Constant-trial value minus the eigenvalue, evaluated directly as the relaxation term. |
| `relaxation_over_rho_squared` | Relaxation term divided by rho squared. |
| `geometric_correction` | Signed correction -((1+epsilon)^(1-2s)-1)/(s(1-2s)). |
| `minus_first_correction_over_scale` | Relaxation minus geometric correction, divided by epsilon for s <= 1/4 or by rho squared for s > 1/4. |
| `predicted_scaled_limit_galerkin` | 1/s below s = 1/4, 1/s + S_galerkin at s = 1/4, and S_galerkin above s = 1/4. |
| `corrector_H1_error` | H1-norm of the finite-space correction discrepancy computed by `run_case`, not a certified continuum error. |
| `schur_residual` | Absolute residual of the scalar eigenvalue Schur identity. |

The three `manuscript/figures/scales_*.dat` files correspond to s = 0.15, 0.25 and 0.35 and have twelve rows each. Their columns are `epsilon`, `geometric`, `spectral`, `combined` and `limit`. The geometric column is minus `geometric_correction` divided by the plotting scale; the spectral column is `trial_minus_lambda` divided by that scale. The combined column is their sum up to rounding, and limit is `predicted_scaled_limit_galerkin`. The scale is epsilon for s <= 1/4 and rho squared otherwise. Text output retains twelve significant digits; the JSON records the full computed precision.

## Graph and kernel results

`verification/results.json` records software versions, the fixed random seed, and graph, tangential-normalization and collar-integral checks. `graph.parameters` contains the fixed positions, weights, coupling/fractional parameters, selected exterior columns, resolvent shift and semigroup time. The two gap fields are graph spectral gaps. The three identity-error fields check elimination, minimized energy and the constant nullvector.

The six entries in `graph.rows` concern perturbations tau B. `tau` is the scale, `eta` the defect bound, `lambda1` the principal eigenvalue, and `first_order_a` the constant-vector first-order approximation. `quadratic_correction_magnitude` is the positive second-order correction magnitude. `second_order_remainder` is lambda1 minus first_order_a plus this magnitude. Fields ending in `_error` and `_error_bound` compare computed errors with explicit finite-graph estimates. The remaining fields record the relative form defect, ground-state bound, remainder divided by tau cubed, maximum spectral shift and split Schur identity error. Precise formulas and matrix norms are in `verification/verify.py`.

The `tangential_normalization` and `collar` groups contain auxiliary integral checks with their own parameter, exact-value, quadrature and error fields. The top-level `status` describes the outcome of program assertions, not certification of the paper's continuum results.
