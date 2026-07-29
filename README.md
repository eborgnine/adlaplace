# adlaplace

Laplace approximations for hierarchical models using [CppAD](https://coin-or.github.io/CppAD/) automatic differentiation (via the **RCppAD** R package) and trust-region optimization.

This repository is a monorepo. The packages most users need are:

| Package | Role |
|---------|------|
| [`RCppAD`](RCppAD/) | Headers-only CppAD for `LinkingTo` / `Imports` (no system CppAD) |
| [`adlaplace`](adlaplace/) | Core Laplace / AD engine |
| [`adlaplaceHgp`](adlaplaceHgp/) | IWP / hierarchical GP terms (`iwp`, `hiwp`, …) |
| [`hpolcc`](hpolcc/) | Hierarchical case-crossover models (`hnlm`, `dirichlet_multinom`) |

Related packages in the same tree: `adlaplaceExample` (custom-density backend example), `adlaplaceFem` (FEM Matérn), `admvn` (MVN / SUN AD utilities). All compiled packages that use CppAD depend on **RCppAD**.

## Install from R-universe (binaries)

Pre-built packages are published at [eborgnine.r-universe.dev](https://eborgnine.r-universe.dev):

```r
install.packages(
  c("RCppAD", "adlaplace", "adlaplaceHgp", "hpolcc"),
  repos = c("https://eborgnine.r-universe.dev", "https://cloud.r-project.org")
)
```

Other packages in this repo (`adlaplaceExample`, `adlaplaceFem`) can be installed the same way. Prefer this over source installs when a binary is available for your platform.

## Quick check

```r
library(adlaplace)   # load adlaplace before hpolcc (CppAD / OpenMP setup)
library(hpolcc)

# Minimal development bundle (no outer optimization)
td <- data.frame(
  count = c(1, 0, 2, 1, 0, 1),
  hum = rnorm(6),
  year = 2002L,
  region = rep(1:2, each = 3),
  date = rep(1:3, 2)
)
fit_dev <- hnlm(
  dirichlet_multinom(count, by = c("year", "region", "date")) ~
    hum + iid(date),
  data = td,
  config = list(num_threads = 1L, num_groups = 2L),
  for_dev = TRUE
)
joint_log_dens(
  fit_dev$ad_pack,
  c(fit_dev$config$opt$init[1], fit_dev$cache$gamma, fit_dev$config$opt$init[-1])
)
```

Load **`adlaplace` before `hpolcc`** (and before `data.table` if you attach it explicitly) so CppAD’s OpenMP setup runs in the core package first.

## Vignettes

HTML vignettes are published by CI to GitHub Pages:

[https://eborgnine.github.io/adlaplace/](https://eborgnine.github.io/adlaplace/)

Examples:

- [adlaplace overview](https://eborgnine.github.io/adlaplace/adlaplace/adlaplace.html)
- [Case-crossover models](https://eborgnine.github.io/adlaplace/adlaplace/casecrossover.html)
- [Dirichlet–multinomial equivalence](https://eborgnine.github.io/adlaplace/adlaplace/dirichlet_multinom.html)
- [Germany BYM example](https://eborgnine.github.io/adlaplace/adlaplace/germany.html)
- [GAMM examples](https://eborgnine.github.io/adlaplace/adlaplace/gamm.html)
- [FEM Matérn methods](https://eborgnine.github.io/adlaplace/adlaplaceFem/fem_matern_methods.html)
- [Loaloa FEM example](https://eborgnine.github.io/adlaplace/adlaplaceFem/loaloa.html)

Case-crossover / `hnlm` usage is documented in the **adlaplace** `casecrossover` vignette (not a separate `hpolcc` HTML vignette). The site refreshes when the `R-CMD-check` workflow completes successfully (Sunday schedule or manual dispatch). Set the repo **Pages** source to **GitHub Actions** once if it is not already.

## System requirements

`adlaplace` and backends need a C++17 compiler. **OpenMP** is recommended (`SystemRequirements: OpenMP`) for multi-threaded fits; install still succeeds without it (`num_threads` is ignored).

CppAD is **not** a system library here. Install the R package [`RCppAD`](RCppAD/) first; packages use `LinkingTo: RCppAD` and `#include <cppad/cppad.hpp>`. There is no need for Homebrew/`apt` CppAD or `libcppad_lib`.

### macOS (Homebrew)

```bash
brew install libomp
```

Xcode command-line tools must be installed (`xcode-select --install`). `configure` detects Homebrew `libomp` when present.

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y libomp-dev g++
```

### Troubleshooting

- **Missing `cppad/cppad.hpp`** — install **RCppAD** first (`R CMD INSTALL RCppAD`), then reinstall packages that `LinkingTo` it. Do not install system `cppad` / `libcppad-dev` for this repo.
- **OpenMP link errors on macOS** — install `libomp` via Homebrew; reinstall `adlaplace` after that so `configure` picks up `-Xpreprocessor -fopenmp` and the `libomp` library. Without OpenMP the package still installs; multi-threading is disabled.
- **`hpolcc` fails to link `adlaplace.so`** — install `adlaplace` into the same R library first, then reinstall `hpolcc`.
- **Verbose configure** — set `ADLAPLACE_CONFIGURE_VERBOSE=1` in the environment before installing.

## License

All packages in this repository are licensed under **GPL (>= 2)**.
See each package’s `DESCRIPTION`. Upstream CppAD is also available under
EPL-2.0 from its own source distribution.
