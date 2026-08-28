# Hierarchical B-splines section for adlaplace.Rmd
#
# Insert after "## Continuous model and FEM weights" (line ~387), before
# "## Range and marginal SD parameterization".
#
# Also update line ~360 to mention optional local refinement.

## Hierarchical B-splines {#sec:hbspline}

The tensor-product construction above uses a single rectangular knot grid.
When data are concentrated in subregions of a padded domain, \pkg{adlaplaceFem}
accepts a **list of refinement levels**---each a \code{SpatRaster} or list of
rasters---via \code{matern(geometry, knots = list(raster0, raster1, \ldots))}.
Level~0 defines a coarse grid over the buffered domain; later levels add knot
lines inside rectangular subregions where finer resolution is needed.

Let $\Omega^0 \supseteq \Omega^1 \supseteq \cdots \supseteq \Omega^L$ denote
the nested refinement regions and $V^\ell$ the tensor-product B-spline space on
the level-$\ell$ global knot grid. A truncated hierarchical basis
\citep{forsey1988hierarchical,kraft1997adaptive,giannelli2012thb} keeps a
level-$\ell$ function $\psi^\ell_{ij}$ active when
$\operatorname{supp}\psi^\ell_{ij}\subseteq\Omega^\ell$ and
$\operatorname{supp}\psi^\ell_{ij}\not\subseteq\Omega^{\ell+1}$. Every active
function embeds in the finest space $V^L$: with sparse $\mathbf{S}$,
$\boldsymbol{\psi}_H=\boldsymbol{\psi}_L\mathbf{S}$, so
\begin{equation}
\mathbf{C}_H = \mathbf{S}^\top \mathbf{C}_L \mathbf{S},\quad
\mathbf{G}_H = \mathbf{S}^\top \mathbf{G}_L \mathbf{S},\quad
\mathbf{G}_{2H} = \mathbf{S}^\top \mathbf{G}_{2L} \mathbf{S},\quad
\mathbf{G}_{3H} = \mathbf{S}^\top \mathbf{G}_{3L} \mathbf{S}.
\end{equation}
The weight precision $\mathbf{Q}_2$ and $\mathbf{Q}_3$ and the \code{fem\_logdet}
atomic are unchanged; only the cached Grams use $\mathbf{C}_H,\ldots$ instead of
a uniform-grid $\mathbf{C}_L,\ldots$. Unlike irregular triangulations in
\pkg{INLA} \citep{lindgren2011spde,rue2009inla}, the hierarchical basis remains
tensor-product on rectangles, so Kronecker assembly applies on the finest level.

**Basis dimension.** A uniform fine grid over the padded domain with $n_x n_y$
B-splines can be replaced by roughly $n_x n_y$ active hierarchical weights
concentrated near data, while retaining coarse resolution elsewhere---for
example the loaloa application (Section~\ref{sec:loaloa}) uses ${\sim}300$ active
weights versus ${\sim}960$ for a uniform $50\,\mathrm{km}$ grid over the same
buffer. Padding by about half a correlation range needs resolution only on the
coarsest level outside the study region.
