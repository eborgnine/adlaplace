# AGENTS.md

This repository contains R packages for Laplace approximations using automatic differentiation.

## Project Structure

- `adlaplace/` — Main R package (Laplace approximations with CppAD)
- `adlaplaceExample/` — Example/extra package (depends on adlaplace)
- `hpol_adlaplace/` — hpolcc package (depends on adlaplace)
- `applications/` — Example applications (canada, india, synthetic)
- `docs/` — Documentation
- `.continue/` — Continue environment configuration

## Build System

Use `make` to build packages:

```bash
make            # Build all packages (adlaplace, adlaplaceExtra, hpolcc)
make adlaplace  # Build only the main package
make clean      # Remove .o and .so files from src folders
```

Build order: `adlaplace` → `adlaplaceExtra`, `hpolcc` (dependencies enforced by Makefile)

## System Requirements

- R (>= 3.6)
- CppAD (automatic differentiation library)
- OpenMP (for parallel computation)
- Rcpp, roxygen2, trustOptim

Install dependencies (Debian/Ubuntu):
```bash
apt-get install r-base r-base-dev libcppad-dev libeigen3-dev r-cran-rcpp r-cran-roxygen2
```

## R Package Dependencies

- Runtime: Matrix, Rcpp, parallel, trustOptim
- Build: Rcpp, trustOptim (LinkingTo)
- Suggests: knitr, rmarkdown, numDeriv, methods, RSpectra, abind

## Development Workflow

1. Edit R or C++ code in package directories
2. Run `make <package>` to compile attributes, generate docs, and build tarball
3. Install with `R CMD INSTALL <package>_<version>.tar.gz`

The Makefile automatically runs:
- `Rcpp::compileAttributes()` — generates RcppExports
- `roxygen2::roxygenize()` — generates NAMESPACE and man pages
- `R CMD build` — creates installable tarball

## Code Style

- C++ code uses Rcpp for R integration
- R documentation via roxygen2 comments
- License: MPL-2.0
