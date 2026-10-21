---
name: mGP near-monotone FEM term
overview: Port Ziang Zhang's near-monotone mGP smoother (FEM method only) into adlaplaceFem as an `mgp` model_term class using adlaplace's existing `random_mult` kernel, with a PSD-calibrated sd prior and a London PM10 case-crossover vignette.
todos:
  - id: fem-core
    content: Port FEM precision/design (one- and two-sided) to R/mgp-fem.R using splineDesign + Gauss-Legendre, no fda dependency
    status: completed
  - id: psd-generic
    content: adlaplace - psd() S4 generic, iwp method (knot-spacing scaled), psd_prior() helper returning a prior() term; bump adlaplace to 0.7.5
    status: completed
  - id: psd-mgp
    content: adlaplaceFem - psd() method for mgp (Green's-function quadrature)
    status: completed
  - id: classes
    content: mgp and mgp_boundary S4 classes + constructor + design/precision/random_info/theta_info/beta_info/replace_vars methods
    status: completed
  - id: tests
    content: Reference tests against Ziang's code (copied helper, skip without fda), random_mult tape checks, small end-to-end fits
    status: completed
  - id: vignette
    content: vignettes/london_mgp.Rmd - case-crossover on adlaplace::london with mgp(pm10) vs iwp(pm10)
    status: completed
  - id: package-meta
    content: DESCRIPTION (Ziang as author, Collate, Description text, version 0.7.2), roxygen, README, R CMD check
    status: completed
---

# mGP near-monotone smoothing in adlaplaceFem (FEM only)

## Sources

- Workflowr repo: <https://github.com/AgueroZZ/near_monotone_smoothing>
  (math in `analysis/index.Rmd`, `analysis/PSD.rmd`; application in
  `analysis/casecross.rmd`; FEM reference in `code/02-FEM.R`).
- The production implementation lives in BayesGP branch
  `bayesgp-near-monotone-usability` (<https://github.com/Bayes-GP/BayesGP>):
  - `R/05_near_monotone_fem.R`: `Compute_Prec()`, `Compute_Prec_rev()`, `Compute_Design()`
  - `R/08_near_monotone_instance_builders.R`: `build_nearmono_instance()`,
    `build_two_sided_mgp_fem()`, `resolve_nearmono_shift()`, domain validation
  - `R/07_near_monotone_helpers.R`: `make_boxcox_transform()`,
    `boundary_basis_matrix()`, `convert_nearmono_sd_prior()`, `resolve_reference_location()`
  - `R/04_near_monotone_state_space.R`: `PSD_compute_mgp()` (via closed-form `SS_cov`),
    `prior_conversion_mgp()`
- Local reference clones: `sandbox/ext/near_monotone_smoothing`, `sandbox/ext/BayesGP`
  (untracked; do not commit).

Out of scope: state-space computation, tIWP2, estimating `a` (BOSS).

## Model

\(L_a f = f'' - \alpha f' = \sigma\xi\) with \(\alpha(x) = -1/(a(x+c))\). The base
model (null space) is \(\{1, m(x)\}\), with \(m\) the Box-Cox transform of \(x+c\),
\(\lambda = (a-1)/a\): \(a=2\) gives \(\sqrt{x+c}\), \(a=1\) gives \(\log(x+c)\),
\(a=\infty\) gives a linear base model (the IWP2).

FEM: \(f = \sum_j w_j\psi_j\) with cubic B-splines, \(w \sim N(0, \sigma^2 Q^{-1})\) and
\(Q_{ij} = \langle L\psi_i, L\psi_j\rangle = G - M + C\), where
\(G = \langle\alpha\psi_i', \alpha\psi_j'\rangle\),
\(M = \langle\alpha\psi_i', \psi_j''\rangle + \text{transpose}\), and
\(C = \langle\psi_i'', \psi_j''\rangle\).

Conventions carried over from BayesGP:

- Reference location `x0` (`"left"` by default, or `"middle"`, `"right"`, or a number).
  Work in the centred coordinate `u = x - x0` with internal shift `cs = x0 + c`, so
  that `c` stays an original-scale shift.
- Default `c`: 0 when the domain minimum is > 0; otherwise `-min + 1e-6 * scale`.
  Require `x + c > 0` over the whole region.
- Positive side: on `[0, max u]`, use `k` equally spaced cubic B-splines (fda convention:
  `k - 2` equally spaced breakpoints) and drop the first two, so `f(x0) = f'(x0) = 0`,
  with \(\alpha(u) = -1/(a(u+cs))\).
- Negative side (only when `x0` is interior or right): `v = x0 - x` with the same
  construction and \(\alpha_{rev}(v) = 1/(a(cs - v))\) (Ziang's `Compute_Prec_rev`).
  Use `B = cbind(B_pos, B_neg)` and `Q = bdiag(Q_pos, Q_neg)`.
- Boundary column (normalized, slope 1 at `x0`, zero at `x0`):
  \(m(u) = ((u+cs)^\lambda - cs^\lambda)/(\lambda cs^{\lambda-1})\), or
  \(cs\log((u+cs)/cs)\) when \(\lambda=0\).
- Because `f(x0) = 0` structurally, the term needs no intercept. This is what the
  case-crossover model (strata absorb the intercept) requires.

## Design decisions

1. **No new C++.** The prior density is \(-\tfrac12\gamma^\top Q\gamma/\sigma^2 - n\log\sigma + \tfrac12\log|Q| - \tfrac n2\log 2\pi\),
   with a fixed Q and one theta. That is exactly adlaplace's `random_mult` kernel
   (`adlaplace/src/ad_pack_random_mult.cpp`), which uses the quadform atomic and
   precomputed `log_det` and `rank`. The `mgp` term sets `density = "random_mult"` and
   `package = "adlaplace"`, so `ad_pack_ptr()` dispatches to
   `adlaplace:::create_ad_shard_random_mult`. `precision()` returns
   `list(Q = <dsCMatrix upper>, log_det, rank = ncol(Q))`; `build_random_shards()`
   already passes list payloads through. `extra_density()` returns `NULL`, because the
   log-determinant is already inside `random_mult`.
   - The FEM kernels (`random_fem_ssq_*`) are not usable: they hard-code two thetas
     (range, sd).
2. **Quadrature.** Replace Ziang's left-point Riemann sum (`accuracy = 0.01` grid) with
   Gauss-Legendre per knot span (`gauss_legendre()`/`gl_on_interval()` in
   `R/grams_1d.R`), using `n_quad = 10` by default. The integrand is rational, not
   polynomial, so `n_quad` stays configurable. Basis evaluation uses `bspline_eval()`
   (splines::splineDesign) with derivs 0/1/2. This removes the fda dependency.
3. **Boundary as a separate term.** Mirror `iwp()` returning `list(iwp, rpoly)`:
   `mgp()` returns `list(<mgp>, <mgp_boundary>)`. `mgp_boundary` copies `rpoly`'s
   pattern: `model_role = "random"`, `sd = Inf` by default (flat, no precision shard), or
   a finite `boundary_sd` (`random_diagonal`, known sd; BayesGP's default prec 0.001 is
   sd ≈ 31.6). Keeping it random means `rmvnldl()`/`sim_random()` include it in the
   simulated curves.
4. **Generic `psd()` in adlaplace.** The predictive SD
   \(\sigma_x(h) = \mathrm{SD}[f(x+h)\mid f(x), f'(x), \dots, f^{(p-1)}(x)]\) is a
   property of the term class, so it becomes an S4 generic next to `design()` and
   `precision()` in `adlaplace/R/classes.R` (or a new `adlaplace/R/psd.R`):
   ```r
   setGeneric("psd", function(term, h, x = NULL, sd = 1, ...) standardGeneric("psd"))
   ```
   - It is vectorised over `h` (and over `x` where `x` matters). It returns `sd` times
     the unit PSD. `sd` defaults to 1 so the same call gives the prior correction
     factor; passing a fitted or simulated sd gives the PSD on the data scale.
   - The default `model_term` method stops with "no psd() method for class ...".
   - **`iwp` method (in adlaplace).** The stored basis is divided by \(s^p\) and the
     precision weights by \(s\), where \(s\) is `iwp_mean_spacing(knots - ref_value)`.
     So the sd theta is on knot-spacing units, and
     \(\mathrm{psd} = \sigma\sqrt{(h/s)^{2p-1}/\big((2p-1)((p-1)!)^2\big)}\).
     That is BayesGP's `PSD_compute_iwp` with `h` replaced by `h/s`, and it is
     invariant to rescaling `x` and the knots together. `x` is ignored because the
     process is stationary in its increments.
   - **`mgp` method (in adlaplaceFem):** see decision 5. `x` is on the original scale
     and defaults to the reference location `x0`. Error if `x + c <= 0` or
     `x + h + c <= 0`.
   - **`psd_prior(term, u, alpha = 0.5, h, x = NULL)` (in adlaplace).** This is
     Ziang's `sd.prior = list(u, alpha, h, x)`: an exponential PC prior with
     \(P(\sigma_x(h) > u) = \alpha\). It returns
     `prior(theta = <sd theta label>, rate = -log(alpha) * psd(term, h, x, sd = 1) / u)`,
     which works with any class that has a `psd()` method. It accepts either a term or
     the list returned by `iwp()` / `mgp()`, and picks the element that has a theta.
   - adlaplaceFem imports `psd` from adlaplace (add it to the `importFrom(adlaplace, ...)`
     block in `000.R`) and calls `exportMethods(psd)`.
5. **mGP PSD by quadrature, not closed forms.**
   \(\mathrm{Var}[f(x+h)\mid f(x), f'(x)] = \int_x^{x+h} \big((g(x+h)-g(\tau))/w(\tau)\big)^2 d\tau\),
   with \(w = (u+cs)^{-1/a}\) and \(g = \int w\). Compute it with `stats::integrate` for
   any `a` (including `a = 1` and `Inf`), instead of porting the six `SS_cov_c*`
   symbolic branches. The integrand is the same on the negative side of a two-sided
   reference, using \(\alpha_{rev}\) with `h` pointing away from `x0`. Tests check it
   against Ziang's `SS_cov`. Port the closed forms only if the quadrature disagrees.
6. `a = Inf` is allowed (\(\alpha = 0\), `Q = C`, `m(u) = u`). This gives an IWP2
   comparison inside the same class.

## Files

New `adlaplaceFem/R/mgp-fem.R` (internal, ported from Ziang with attribution in roxygen):

- `mgp_shift(c, domain)`, `mgp_reference(x, ref_value)`, `mgp_check_domain()`
- `mgp_boxcox(u, a, cs)`: the normalized boundary transform
- `mgp_side_basis(region_max, k)`: open knot vector and dropped columns
- `mgp_side_precision(knots, a, cs, side, n_quad)`: `G - M + C` via GL quadrature
- `mgp_side_design(u, knots, side)`
- `mgp_fem(domain, x0, a, cs, k, n_quad)`: assembles both sides and returns knots, `Q`,
  `log_det` (sparse Cholesky), and `rank`

New `adlaplaceFem/R/mgp.R`:

- `setClass("mgp", contains = "model_term")` with slots `a`, `c`, `shift`, `region`,
  `k`, `n_quad`, `fem` (cached list), `two_sided`. Prototype: `density = "random_mult"`,
  `ad_kind = "random"`, `package = "adlaplace"`. Set `knots = region breakpoints` so
  `sim_random()` builds its grid over the region.
- Constructor `mgp(x, a = 2, c = NULL, k = 30, region = NULL, ref_value = "left",
  boundary_sd = Inf, n_quad = 10, init, lower, upper, parscale, log = TRUE)`.
  `region` is required for now (the term is built before data is seen; this matches
  `iwp()` needing knots). Validate `x + c > 0` on the region.
- Methods: `design` (sparse B at `data[[x]]`, error outside region, colnames
  `<x>_mgp_k##`), `precision` (random_mult list), `random_info`, `theta_info` (one
  row, `<x>_mgp` sd), `beta_info` (NULL), `replace_vars`.
- `setClass("mgp_boundary")` + `design` / `precision` / `random_info`, mirroring `rpoly`.
- `setMethod("psd", "mgp", ...)`, using the internal `mgp_psd_unit(x, h, a, cs, side)`.

Update `NAMESPACE` (roxygen), `Collate`, and `000.R` imports (`psd`).

New `adlaplace/R/psd.R`:

- The `psd` generic, its `model_term` default (error), and the `iwp` method.
- `psd_prior()` (exported). Add to `Collate` / NAMESPACE, and document under
  `?psd` with a short section in the `iwp` docs.
- Bump adlaplace to 0.7.5; adlaplaceFem then `Depends: adlaplace (>= 0.7.5)`.

## Tests (`tests/testthat/test-mgp.R`, helper `helper-mgp-reference.R`)

The helper is a verbatim copy of Ziang's `Compute_Prec`, `Compute_Prec_rev`,
`Compute_Design`, `boundary_basis_matrix`, `make_boxcox_transform`, and `SS_cov*`, with
a source comment. FEM comparisons use `skip_if_not_installed("fda")`; fda goes in
Suggests.

- Q (one- and two-sided, `a ∈ {2, 1, -1, 3}`) matches Ziang's `Compute_Prec` at
  `accuracy = 1e-4`, to a relative tolerance of about 1e-3. The design matrix matches
  `Compute_Design` exactly.
- The boundary column matches `boundary_basis_matrix(normalized = TRUE)` for
  left/middle/right reference locations.
- `psd(<mgp>, h, x)` matches `sqrt(SS_cov(x + h, x, c, a))` for `a ∈ {2, 1, -1, 3, -3}`.
- Self-consistency for `mgp`: because `f(x0) = f'(x0) = 0`,
  `sqrt(diag(B Q^{-1} B'))` at `x0 + h` approaches `psd(term, h, x = x0)` as `k` grows.
- `psd_prior()` on an `mgp` term gives the same `rate` as Ziang's
  `prior_conversion_mgp()`.

In `adlaplace/tests/testthat/test-psd.R`:

- `psd(<iwp>)` equals BayesGP's `PSD_compute_iwp(h, p)` when the mean knot spacing is 1,
  for `p ∈ {1, 2, 3}`.
- Scale invariance: multiplying `x`, the knots and `h` by 10 leaves `psd()` unchanged.
- Self-consistency for `iwp`: with `x = ref_value` and `ref + h` a knot,
  `sd^2 * diag(B Q^{-1} B')` from the term's own `design()`/`precision()` equals
  `psd()^2` exactly. The O-spline is exact at knots, and the polynomial part is the
  separate `rpoly` term.
- `psd_prior()` resolves the right theta label for `iwp()` lists and errors on classes
  without a `psd()` method.
- FEM prior covariance `B Q^{-1} B'` at a few points approaches Ziang's
  `mspline_cov` as `k` grows (loose tolerance, `k = 80`).
- Tape: build a `random_mult` shard from `mgp`. `joint_log_dens` and its gradient
  equal the dense R density, and the gradient matches numDeriv.
- End-to-end: (a) `normal(y) ~ mgp(x, a = 2) + ...` on data simulated from
  \(\sqrt{x}\) plus noise recovers the curve; (b) a london-2003 case-crossover fit runs
  and returns finite `logLik`.
- `replace_vars()` and `sim_random()` work with an `mgp` term.

## Vignette `adlaplaceFem/vignettes/london_mgp.Rmd`

Follows `adlaplace/vignettes/casecrossover.Rmd`, including its `vignette_full` switch
(the year-2003 subset for fast builds).

1. Background: near-monotone exposure–response; the base model \(\sqrt{x+c}\); mGP as
   deviation from that base; cite Ziang's repo and Ramsay (1998).
2. Data: `data("london", package = "adlaplace")`, dropping the 1 NA `pm10`. PM10
   ranges from 8.1 to 99.6.
3. Model:
   ```r
   dirichlet_multinom(all, by = c("year", "month", "dow"), init = 0.02) ~
     0 + rh +
     mgp(pm10, a = 2, c = 0, k = 30, region = c(8, 100), ref_value = "left") +
     iwp(tmean, p = 3, knots = seq(-2, 30, by = 2), ref_value = 20) +
     psd_prior(pm10_mgp, u = 0.01, alpha = 0.5, h = 5, x = 8)
   ```
   This is Ziang's PSD prior from `casecross.rmd`; `pm10_mgp` is the same `mgp(...)`
   call built beforehand. Check that the adlaplace formula parser evaluates a
   pre-built term object; otherwise build the term inline twice. The `iwp(pm10)`
   comparison fit uses `psd_prior(<iwp>, u = 0.01, alpha = 0.5, h = 5)`, so both
   smoothers get matched priors on the 5-unit PSD.
4. Fit with `adlaplace()`; print `par_info`.
5. Plot relative risk vs PM10 (`rmvnldl` + `sim_random`, `transform = exp`) for the mGP
   fit, overlaid with the `iwp(pm10, p = 2)` fit from the adlaplace vignette re-based
   to the same reference.
6. Sensitivity: `a ∈ {1, 2, Inf}` with a `logLik` comparison table, plus a two-sided
   fit with `ref_value = 40` (eval only when `vignette_full`).
7. Prior vs posterior of \(\sigma_x(5)\) on the PSD scale for both smoothers, using
   `psd(term, h = 5, x = ..., sd = <fitted sd>)`. Also plot the mGP PSD against `x`,
   which shows that the mGP is not stationary while the IWP is.

## Package metadata

- DESCRIPTION: add
  `person("Ziang", "Zhang", email = "ziangzhang@uchicago.edu", role = c("aut"))`,
  with a comment noting the near-monotone mGP FEM. Extend Description; add `mgp.R` and
  `mgp-fem.R` to Collate; add `fda` and `numDeriv` to Suggests (tests only); bump to
  0.7.2.
- README: one line on `mgp()`, and one on `psd()` / `psd_prior()` in adlaplace's README.
- adlaplace: add `psd.R` to Collate, regenerate NAMESPACE, add a NEWS entry, and run its
  test suite before rebuilding adlaplaceFem against it.
- `devtools::document()`, `devtools::test()`, `R CMD build` + `R CMD check --no-manual`.

## Risks / checks during implementation

- `random_mult` has so far only been exercised through the low-level germany pipeline,
  not through formula `model_data()`. Check gamma-map compaction, the
  `theta_map_for_term` grep on `<x>_mgp` (it must not also match `<x>_mgp_boundary`),
  and the Hessian-inverse trace path (`hinv_trace_index`) for analytic outer gradients.
- Conditioning: the dropped-column cubic basis makes Q PD, but with `cs` small
  (`c = 0`, left reference at the region minimum) the \(\alpha^2\) term is large near
  `u = 0`. Check `log_det` stability and possibly default `n_quad` higher on the first
  span.
- Basis convention: Ziang's `region` for each side is `[0, max |u|]` from the data plus
  the region. Here the term owns `region`, so check that `ref_value = 40` with
  `region = c(8, 100)` gives sides of length 32 and 60 with `k` basis functions each,
  as in BayesGP.
