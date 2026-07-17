# adlaplace

Laplace approximations for hierarchical models using [CppAD](https://coin-or.github.io/CppAD/) automatic differentiation (via the **RCppAD** R package) and trust-region optimization.

This repository is a monorepo. The packages most users need are:

| Package | Role |
|---------|------|
| [`RCppAD`](RCppAD/) | Headers-only CppAD for `LinkingTo` / `Imports` (no system CppAD) |
| [`adlaplace`](adlaplace/) | Core Laplace / AD engine |
| [`adlaplaceHgp`](adlaplaceHgp/) | IWP / hierarchical GP terms (`iwp`, `hiwp`, …) |
| [`hpolcc`](hpolcc/) | Hierarchical case-crossover models (`hnlm`, `dirichlet_multinom`) |

Related packages in the same tree: `adlaplaceExample` (custom-density backend example), `adlaplaceGrf` (FEM Matérn), `admvn` (MVN / SUN AD utilities). All compiled packages that use CppAD depend on **RCppAD**.

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

### R packages

Install CRAN dependencies first (or let `devtools` / `remotes` pull them):

```r
install.packages(c(
  "Matrix", "Rcpp", "RcppEigen", "trustOptim",
  "data.table", "timeDate"
))
```

Then install **RCppAD** from this repo (not yet assumed on CRAN in local workflows).

## Install from a local clone

Clone the repo, then install **in this order** (`RCppAD` before anything that compiles against CppAD):

```bash
git clone https://github.com/eborgnine/adlaplace.git
cd adlaplace
```

```r
# From the repository root
devtools::install("RCppAD", upgrade = "never")
devtools::install("adlaplace", upgrade = "never")
devtools::install("adlaplaceHgp", upgrade = "never")
devtools::install("hpolcc", upgrade = "never")
```

Or with `R CMD INSTALL`:

```bash
R CMD INSTALL RCppAD
R CMD INSTALL adlaplace
R CMD INSTALL adlaplaceHgp
R CMD INSTALL hpolcc
```

`configure` runs automatically during install and writes `src/Makevars` from `src/Makevars.in` (OpenMP detection only; CppAD comes from RCppAD).

## Install from GitHub

Using [`remotes`](https://remotes.r-lib.org/) (or `devtools`):

```r
install.packages("remotes")

remotes::install_github("eborgnine/adlaplace", subdir = "RCppAD")
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplace")
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplaceHgp")
remotes::install_github("eborgnine/adlaplace", subdir = "hpolcc")
```

Install `RCppAD` then `adlaplace` before the others so LinkingTo / runtime linking can find headers and `adlaplace.so`.

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
  config = list(num_threads = 1L, num_shards = 2L),
  for_dev = TRUE
)
joint_log_dens(
  fit_dev$ad_fun,
  c(fit_dev$config$opt$init[1], fit_dev$cache$gamma, fit_dev$config$opt$init[-1])
)
```

Load **`adlaplace` before `hpolcc`** (and before `data.table` if you attach it explicitly) so CppAD’s OpenMP setup runs in the core package first.

## Optional packages

These also need **RCppAD** (and usually **adlaplace**) already installed:

```r
# Custom observation densities (skew-normal backend demo)
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplaceExample")

# FEM Matérn spatial fields
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplaceGrf")

# MVN / SUN AD utilities
remotes::install_github("eborgnine/adlaplace", subdir = "admvn")
```

## Refreshing CppAD (maintainers)

Upstream headers are vendored under [`RCppAD/inst/include/cppad/`](RCppAD/inst/include/cppad/). To update:

```bash
cd RCppAD
./tools/update-cppad.sh 20260000.0   # or a path to include/
R CMD INSTALL .
```

See [`RCppAD/README.md`](RCppAD/README.md). Do **not** download CppAD at `R CMD INSTALL` time.

## Troubleshooting

- **Missing `cppad/cppad.hpp`** — install **RCppAD** first (`R CMD INSTALL RCppAD`), then reinstall packages that `LinkingTo` it. Do not install system `cppad` / `libcppad-dev` for this repo.
- **OpenMP link errors on macOS** — install `libomp` via Homebrew; reinstall `adlaplace` after that so `configure` picks up `-Xpreprocessor -fopenmp` and the `libomp` library. Without OpenMP the package still installs; multi-threading is disabled.
- **`hpolcc` fails to link `adlaplace.so`** — install `adlaplace` into the same R library first, then reinstall `hpolcc`.
- **Verbose configure** — set `ADLAPLACE_CONFIGURE_VERBOSE=1` in the environment before installing.

## License

See each package’s `DESCRIPTION`:

- `RCppAD` — EPL-2.0 or GPL-2+ (`License: file LICENSE`; vendored CppAD)
- `adlaplace` — MPL-2.0
- `hpolcc`, `adlaplaceHgp`, and most backends — GPL-3
