PKGS := adlaplace adlaplaceExample hpolcc
ADLAPLACE_DIR := adlaplace
ADLAPLACE_EXAMPLE_DIR := adlaplaceExample
HPOLCC_DIR := hpolcc

.DEFAULT_GOAL := all

.PHONY: all clean $(PKGS) dirichlet_multinom.pdf

all: $(PKGS) dirichlet_multinom.pdf

clean:
	@echo "==> Cleaning .o and .so files under src folders"
	@find . -type f \( -name '*.o' -o -name '*.so' \) -path '*/src/*' -delete



adlaplace:
	@echo "==> Running cleanup for adlaplace"
	cd $(ADLAPLACE_DIR) && ./cleanup
	@echo "==> Running compileAttributes for adlaplace from $(ADLAPLACE_DIR)"
	Rscript -e "Rcpp::compileAttributes('$(ADLAPLACE_DIR)')"
	@echo "==> Running roxygen2 for adlaplace from $(ADLAPLACE_DIR)"
	Rscript -e "roxygen2::roxygenize('$(ADLAPLACE_DIR)')"
	@echo "==> Building package adlaplace from $(ADLAPLACE_DIR)"
	R CMD build --no-build-vignettes $(ADLAPLACE_DIR)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(ADLAPLACE_DIR)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(ADLAPLACE_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

adlaplaceExample: adlaplace
	@echo "==> Running cleanup for adlaplaceExample"
	cd $(ADLAPLACE_EXAMPLE_DIR) && ./cleanup
	@echo "==> Running compileAttributes for adlaplaceExample from $(ADLAPLACE_EXAMPLE_DIR)"
	Rscript -e "Rcpp::compileAttributes('$(ADLAPLACE_EXAMPLE_DIR)')"
	@echo "==> Running roxygen2 for adlaplaceExample from $(ADLAPLACE_EXAMPLE_DIR)"
	Rscript -e "roxygen2::roxygenize('$(ADLAPLACE_EXAMPLE_DIR)')"
	@echo "==> Building package adlaplaceExample from $(ADLAPLACE_EXAMPLE_DIR)"
	R CMD build --no-build-vignettes $(ADLAPLACE_EXAMPLE_DIR)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(ADLAPLACE_EXAMPLE_DIR)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(ADLAPLACE_EXAMPLE_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

hpolcc: adlaplace
	@echo "==> Running cleanup for hpolcc"
	cd $(HPOLCC_DIR) && ./cleanup
	@echo "==> Skipping roxygen2 for hpolcc (requires adaplace dependency)"
	@echo "==> Building package hpolcc from $(HPOLCC_DIR)"
	R CMD build --no-build-vignettes $(HPOLCC_DIR)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

dirichlet_multinom.pdf: hpolcc/vignettes/dirichlet_multinom.Rmd
	pandoc hpolcc/vignettes/dirichlet_multinom.Rmd -o dirichlet_multinom.pdf