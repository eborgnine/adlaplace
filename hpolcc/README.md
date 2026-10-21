# hpolcc

Hierarchical functional case-crossover models for air-pollution health effects,
built on [`adlaplace`](https://github.com/eborgnine/adlaplace) and
[`adlaplaceHgp`](https://github.com/eborgnine/adlaplace).

## Purpose

Fits stratified count models (Dirichlet–multinomial observation layer) with
case-crossover time structure, random date effects, and optional IWP smooths
for exposure–response curves.

## Main entry points

| Function | Role |
|----------|------|
| `hnlm()` | Fit a hierarchical Laplace model (flat `hnlm` object or `for_dev` bundle) |
| `dirichlet_multinom()` | Response family for overdispersed multinomial counts |
| `case_crossover()` | Build case-crossover strata from calendar variables |
| `remove_holidays()` | Filter observation dates to school-year windows and holidays |

Install after `adlaplace` and `adlaplaceHgp`:

```r
devtools::install("hpolcc")
```

See the monorepo [`README.md`](../README.md) and the
[case-crossover vignette](https://eborgnine.github.io/adlaplace/adlaplace/casecrossover.html)
in **adlaplace** for a worked example.
