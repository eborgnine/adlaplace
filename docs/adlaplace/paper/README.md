# Paper additions for hierarchical B-splines

This directory holds material for the JSS paper at
`~/research/healthcanada/docs/adlaplace/paper`. Merge into the main repository as follows:

1. **Section text:** copy [`sections/hierarchical-bsplines.md`](sections/hierarchical-bsplines.md) into `adlaplace.Rmd` after the "Continuous model and FEM weights" subsection (before "Range and marginal SD parameterization"). Update the FEM section opener to mention optional local refinement.

2. **Bibliography:** append entries from [`adlaplace.bib`](adlaplace.bib).

3. **Loaloa example:** replace `replication/04-loaloa.R` and update §7.5 (`loaloa-formula`, `loaloa-run`, add `hbspline-figure` chunk for `plot_loaloa_knots()`).

4. **Covariance validation:** incorporate [`replication/fem-validation-hierarchical.R`](replication/fem-validation-hierarchical.R) into the `fem-run` chunk; add a hierarchical row to Table~\ref{tab:fem} and curve to Figure~\ref{fig:fem-corr}.

5. Rebuild: `make clean && make` in the paper directory (clear `cache/fem-*` and `cache/loaloa-*` first).
