PKGS := adlaplace adlaplaceExample hpolcc admvn
ADLAPLACE_DIR := adlaplace
ADLAPLACE_EXAMPLE_DIR := adlaplaceExample
HPOLCC_DIR := hpolcc

.DEFAULT_GOAL := all

.PHONY: all clean $(PKGS) dirichlet_multinom.pdf

all:
	 R -e "devtools::document(\"adlaplace\")" ;
	 R CMD build --no-build-vignettes adlaplace; 
	 R -e "devtools::document(\"hpolcc\")";
	 R CMD build --no-build-vignettes hpolcc 


clean:
	@echo "==> Cleaning package build artifacts (cleanup scripts + src objects)"
	@for pkg in $(PKGS); do \
		if [ -f "$$pkg/cleanup" ]; then \
			echo "==> $$pkg/cleanup"; \
			(cd "$$pkg" && sh ./cleanup); \
		fi; \
	done
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
	echo "==> Built $$TARBALL (installation skipped)"

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
	echo "==> Built $$TARBALL (skipping installation)"

hpolcc: adlaplace
	@echo "==> Running cleanup for hpolcc"
	cd $(HPOLCC_DIR) && ./cleanup
	@echo "==> Running compileAttributes for hpolcc from $(HPOLCC_DIR)"
	Rscript -e "Rcpp::compileAttributes('$(HPOLCC_DIR)')"
	@echo "==> Running roxygen2 for hpolcc from $(HPOLCC_DIR)"
	Rscript -e "roxygen2::roxygenize('$(HPOLCC_DIR)')"
	@echo "==> Building package hpolcc from $(HPOLCC_DIR)"
	R CMD build --no-build-vignettes $(HPOLCC_DIR)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

admvn:
	@echo "==> Running cleanup for admvn"
	cd admvn && ./cleanup
	@echo "==> Running configure for admvn"
	cd admvn && ./configure
	@echo "==> Running compileAttributes for admvn from admvn"
	Rscript -e "Rcpp::compileAttributes('admvn')"
	@echo "==> Running roxygen2 for admvn from admvn"
	Rscript -e "roxygen2::roxygenize('admvn')"
	@echo "==> Building package admvn from admvn"
	R CMD build --no-build-vignettes admvn
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' admvn/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' admvn/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL (installation skipped)"

dirichlet_multinom.pdf: hpolcc/vignettes/dirichlet_multinom.Rmd
	pandoc hpolcc/vignettes/dirichlet_multinom.Rmd -o dirichlet_multinom.pdf