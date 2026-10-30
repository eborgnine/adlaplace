# Hierarchical term interface changes

Install from the adlaplace repo root with `make` (or at least `make adlaplace adlaplaceHgp adlaplaceFem`). Required versions:

- adlaplace >= 0.7.6
- adlaplaceFem >= 0.7.4
- adlaplaceHgp >= 0.6.2

`hpolcc` is deprecated and was not updated.

## `hmgp()` (new, adlaplaceFem)

Grouped near-monotone mGP. One shared deviation SD across levels of `by`, plus an iid boundary slope per level with its own SD.

```r
adlaplaceFem::hmgp(
  x, a = 2, knots, ref_value = "left", by, levels = NULL, n_quad = 10,
  init = 0.02, lower = 1e-9, upper = Inf, parscale = 1, log = TRUE,
  slope_init = 0.02, slope_lower = 1e-9, slope_upper = Inf,
  slope_parscale = 1, slope_log = TRUE,
  include_global = TRUE, boundary_sd = Inf
)
```

`include_global = TRUE` is the default, so one call returns four terms, in order:

1. `<x>_mgp` — common curve, theta label `<x>_mgp`
2. `<x>_mgp_boundary` — flat global slope (`boundary_sd = Inf` means no theta)
3. `<x>_hmgp` — per-group deviation, theta label `<x>_hmgp`
4. `<x>_hslope_mgp` — per-group slope, theta label `<x>_hslope_mgp`

Drop a separate `mgp()` call, or pass `include_global = FALSE` if the common curve is built on its own (for example with a different `init`). The global term reuses `init` / `lower` / `upper` / `parscale` / `log`.

Gamma labels:

- deviation: `<x>_hmgp_k##_g<level>` (basis fastest within each level)
- slope: `<x>_hslope_mgp_g<level>`

Pass `levels` in the fitted order. A prediction grid then gets the same column order as `sample_gamma()` row names. Missing values in `by`, and levels that were not supplied, are errors.

`psd()` on the `hmgp` term matches `psd()` on `mgp`. `psd_prior()` on the whole list hits the first theta, which is `<x>_mgp`. Pass the `<x>_hmgp` element to put a PC prior on the deviation SD.

## `hiwp()` (breaking, adlaplaceHgp)

`hiwp` now extends `iwp`. Positional vectors of length `2p + 4` are gone.

```r
adlaplaceHgp::hiwp(
  x, p = 2, ref_value = 0, knots, by, levels = NULL,
  init = 0.02, lower = 1e-9, upper = Inf, parscale = 1, log = TRUE,
  slope_init = 0.02, slope_lower = 1e-9, slope_upper = Inf,
  slope_parscale = 1, slope_log = TRUE,
  boundary_is_random = TRUE, include_poly = TRUE, include_global = TRUE
)
```

- `init`, `lower`, `upper`, `parscale`, and `log` are scalars. They set the deviation SD and, when `include_global` is `TRUE`, the global `iwp` SD as well.
- `slope_*` is recycled to `p - 1` and sets each `hrpoly` term.
- A different global SD is a separate `iwp()` plus `hiwp(..., include_global = FALSE)`.
- Default return order is now `<x>_iwp`, `<x>_rpoly` or `<x>_fpoly`, `<x>_hiwp`, then `<x>_hrpoly_<d>`.
- `levels` works as for `hmgp`. Unknown levels and `NA` in `by` are errors.
- Gamma labels are unchanged in form: `<x>_hiwp_k##_g<level>`, `<x>_hrpoly_<d>_g<level>`.

## `rsmatern()` (adlaplaceFem)

Coefficient order is unchanged. Gamma labels changed from `<label>_<level>_b##` to `<label>_b##_g<level>`.

## Exported helpers (adlaplace)

`add_by_levels()`, `group_index()`, `scatter_by_group()`, and `group_gamma_labels()` are exported for other grouped terms. `add_by_levels()` errors on `NA` in the grouping column, and errors (it no longer warns and extends) when data has a level that was not supplied.
