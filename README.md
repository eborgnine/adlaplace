# adlaplace

Laplace approximations for hierarchical models using [CppAD](https://coin-or.github.io/CppAD/) automatic differentiation and trust-region optimization.

This repository is a monorepo. The packages most users need are:

| Package | Role |
|---------|------|
| [`RCppAD`](RCppAD/) | Vendored CppAD headers (`LinkingTo`) |
| [`adlaplace`](adlaplace/) | Core Laplace / AD engine |
| [`adlaplaceHgp`](adlaplaceHgp/) | IWP / hierarchical GP terms (`iwp`, `hiwp`, …) |
| [`hpolcc`](hpolcc/) | Hierarchical case-crossover models (`hnlm`, `dirichlet_multinom`) |

Related packages in the same tree: `adlaplaceExample` (custom-density backend example), `adlaplaceGrf` (FEM Matérn), `admvn` (MVN / SUN AD utilities).

## System requirements

`adlaplace` and backends need a C++17 compiler. **OpenMP** is recommended (`SystemRequirements: OpenMP`) for multi-threaded fits; install still succeeds without it (`num_threads` is ignored). CppAD comes from the R package **RCppAD** (no system CppAD / `libcppad_lib`).

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

## Install from a local clone

Clone the repo, then install **in this order** (`RCppAD` before `adlaplace`; `hpolcc` depends on `adlaplace` and `adlaplaceHgp`):

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

`configure` runs automatically during install and writes `src/Makevars` from `src/Makevars.in`.

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

```r
# Custom observation densities (skew-normal backend demo)
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplaceExample")

# FEM Matérn spatial fields
remotes::install_github("eborgnine/adlaplace", subdir = "adlaplaceGrf")
```

## Troubleshooting

- **Missing `cppad/cppad.hpp`** — install **RCppAD** first (`R CMD INSTALL RCppAD`), then reinstall packages that `LinkingTo` it.
- **OpenMP link errors on macOS** — install `libomp` via Homebrew; reinstall `adlaplace` after that so `configure` picks up `-Xpreprocessor -fopenmp` and the `libomp` library. Without OpenMP the package still installs; multi-threading is disabled.
- **`hpolcc` fails to link `adlaplace.so`** — install `adlaplace` into the same R library first, then reinstall `hpolcc`.
- **Verbose configure** — set `ADLAPLACE_CONFIGURE_VERBOSE=1` in the environment before installing.

## License

See each package’s `DESCRIPTION` (`adlaplace` is MPL-2.0; `hpolcc` and `adlaplaceHgp` are GPL-3).
