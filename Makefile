PKGS := adlaplace adlaplaceExtra hpolcc
ADLAPLACE_DIR := adlaplace
ADLAPLACE_EXTRA_DIR := adlaplaceExample
HPOLCC_DIR := hpolcc

.DEFAULT_GOAL := all

.PHONY: all clean $(PKGS) dirichlet_multinom.pdf

all: $(PKGS) dirichlet_multinom.pdf

clean:
	@echo "==> Cleaning .o and .so files under src folders"
	@find . -type f \( -name '*.o' -o -name '*.so' \) -path '*/src/*' -delete

define build_pkg
	@echo "==> Running compileAttributes for $(1) from $(2)"
	Rscript -e "Rcpp::compileAttributes('$(2)')"
	@echo "==> Running roxygen2 for $(1) from $(2)"
	Rscript -e "roxygen2::roxygenize('$(2)')"
	@echo "==> Generating NAMESPACE from existing exports for $(1)"
	Rscript -e "exports <- c('f', 'hnlm', 'ccDesign', 'cond_sim', 'cond_sim_iwp', 'formatHpolData', 'format_parameters', 'getAdFun_r', 'get_effect', 'merge_data_cc', 'removeUnusedStrata'); cat(paste0('export(', exports, ')\n'), file = '$(2)/NAMESPACE')"
	@echo "==> Building package $(1) from $(2)"
	R CMD build --no-build-vignettes $(2)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(2)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(2)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"
endef

adlaplace:
	$(call build_pkg,adlaplace,$(ADLAPLACE_DIR))

adlaplaceExtra: adlaplace
	$(call build_pkg,adlaplaceExtra,$(ADLAPLACE_EXTRA_DIR))

hpolcc: adlaplace
	@echo "==> Skipping compileAttributes for hpolcc (not a standard R package)"
	@echo "==> Skipping roxygen2 for hpolcc (not a standard R package)"
	@echo "==> Generating NAMESPACE from existing exports for hpolcc"
	Rscript -e "exports <- c('f', 'hnlm'); cat(paste0('export(', exports, ')\n'), file = '$(HPOLCC_DIR)/NAMESPACE')"
	@echo "==> Building package hpolcc from $(HPOLCC_DIR)"
	R CMD build --no-build-vignettes $(HPOLCC_DIR)
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' $(HPOLCC_DIR)/DESCRIPTION | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> Built $$TARBALL"

dirichlet_multinom.pdf: hpolcc/vignettes/dirichlet_multinom.Rmd
	pandoc hpolcc/vignettes/dirichlet_multinom.Rmd -o dirichlet_multinom.pdf
