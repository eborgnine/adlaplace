# adlaplaceHgp

Hierarchical and random-slope Gaussian process terms for use with
[`adlaplace`](https://github.com/eborgnine/adlaplace).

## Purpose

This package extends the `adlaplace` formula / S4 term system with integrated
Wiener process (IWP) smooths and related random-effect structures for
hierarchical models (group-specific smooths, scaled random slopes, etc.).

## Main entry points

| Function | Role |
|----------|------|
| `hiwp()` | Hierarchical IWP smooth with grouping (`by`) |
| `hrpoly()` | Hierarchical random polynomial by group |
| `rsiwp()`, `rsiid()` | Random-slope IWP / intercept terms |
| `rsfpoly()`, `rsrpoly()` | Fixed / random polynomial interaction terms |
| `cond_sim_iwp()`, `cond_sim_iwp_at()` | Conditional simulation of IWP/HIWP curves |

Install after `adlaplace`:

```r
devtools::install("adlaplaceHgp")
```

See the monorepo root [`README.md`](../README.md) for the full install order.
