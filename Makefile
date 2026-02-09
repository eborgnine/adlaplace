PKG_DIR := adlaplace
VERSION := $(shell sed -n 's/^Version:[[:space:]]*//p' $(PKG_DIR)/DESCRIPTION | head -n 1)
TARBALL := $(PKG_DIR)_$(VERSION).tar.gz

.DEFAULT_GOAL := $(TARBALL)

.PHONY: $(TARBALL)

$(TARBALL):
	@echo "==> Running compileAttributes for $(PKG_DIR)"
	Rscript -e "Rcpp::compileAttributes('$(PKG_DIR)')"
	@echo "==> Running roxygenize for $(PKG_DIR)"
	Rscript -e "roxygen2::roxygenize('$(PKG_DIR)', load = 'source')"
	@echo "==> Building package $(PKG_DIR)"
	R CMD build --no-build-vignettes $(PKG_DIR)
	@test -f "$(TARBALL)" || { echo "Expected tarball $(TARBALL) not found"; exit 1; }
	@echo "==> Built $(TARBALL)"
