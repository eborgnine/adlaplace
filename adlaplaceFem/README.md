# adlaplaceFem

Finite-element approximations for hierarchical models with
[adlaplace](../adlaplace). Currently focuses on tensor-product B-spline FEM
ingredients for Matérn Gaussian random fields; additional FEM approximations
will live here.

## Features

- `fem_bspline()` — design matrix `A` and Grams `C`, `G`, `G2` (and `G3` for degree ≥ 3)
- `fem_precision()` / `fem_precision_payload()` — assemble `Q₂`/`Q₃` and chol pattern for `random_fem_ssq_*` / `random_fem_det_*`
- `matern()` — formula term for 2D Matérn GRFs via FEM
- Optional **terra** `SpatRaster` / `SpatVector` coords (Suggests only)

## Vignettes

- `fem_matern_methods` — SPDE/FEM mathematics and Matérn vs FEM covariance comparison
- `fem_joint_density` — joint density / gradient / Hessian via `adlaplace` (no optimization)
- `loaloa` — prevalence example with FEM Matérn
