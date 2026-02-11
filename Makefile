PKGS := hpolcc
HPOLCC_DIR := hpol_adlaplace

.DEFAULT_GOAL := all

.PHONY: all $(PKGS)

all: hpolcc

hpolcc:
	@echo "==> Running compileAttributes for hpolcc"
	Rscript -e "Rcpp::compileAttributes('$(HPOLCC_DIR)')"
	@echo "==> Running roxygenize for hpolcc"
	Rscript -e "roxygen2::roxygenize('$(HPOLCC_DIR)', load = 'source')"
	@echo "==> Building package hpolcc from $(HPOLCC_DIR)"
	R CMD build --no-build-vignettes $(HPOLCC_DIR)
	@VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="hpolcc_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"
