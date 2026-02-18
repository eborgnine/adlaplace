PKGS := adlaplace adlaplaceExtra hpolcc
ADLAPLACE_DIR := adlaplace
ADLAPLACE_EXTRA_DIR := adlaplaceExample
HPOLCC_DIR := hpol_adlaplace

.DEFAULT_GOAL := all

.PHONY: all clean $(PKGS)

all: $(PKGS)

clean:
	@echo "==> Cleaning .o and .so files under src folders"
	@find . -type f \( -name '*.o' -o -name '*.so' \) -path '*/src/*' -delete

define build_pkg
	@echo "==> Running compileAttributes for $(1) from $(2)"
	Rscript -e "Rcpp::compileAttributes('$(2)')"
	@echo "==> Running roxygenize for $(1)"
	Rscript -e "roxygen2::roxygenize('$(2)', load = 'source')"
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
	$(call build_pkg,hpolcc,$(HPOLCC_DIR))
