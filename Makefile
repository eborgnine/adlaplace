PKGS := adlaplace adlaplaceExample

.DEFAULT_GOAL := all

.PHONY: all $(PKGS)

all: adlaplace adlaplaceExample

adlaplace:
	@echo "==> Running compileAttributes for adlaplace"
	Rscript -e "Rcpp::compileAttributes('adlaplace')"
	@echo "==> Running roxygenize for adlaplace"
	Rscript -e "roxygen2::roxygenize('adlaplace', load = 'source')"
	@echo "==> Building package adlaplace"
	R CMD build --no-build-vignettes adlaplace
	@VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' adlaplace/DESCRIPTION | head -n 1)"; \
	TARBALL="adlaplace_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

adlaplaceExample: adlaplace
	@echo "==> Running compileAttributes for adlaplaceExample"
	Rscript -e "Rcpp::compileAttributes('adlaplaceExample')"
	@echo "==> Running roxygenize for adlaplaceExample"
	Rscript -e "roxygen2::roxygenize('adlaplaceExample', load = 'source')"
	@echo "==> Building package adlaplaceExample"
	R CMD build --no-build-vignettes adlaplaceExample
	@VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' adlaplaceExample/DESCRIPTION | head -n 1)"; \
	TARBALL="adlaplaceExample_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"
