
.PHONY: FORCE
FORCE:

# discover package dirs (subdirectories containing DESCRIPTION)
PKGS := $(patsubst %/DESCRIPTION,%,$(wildcard */DESCRIPTION))

.PHONY: $(PKGS)

$(PKGS): FORCE
	@echo "==> Running compileAttributes for $@"
	Rscript -e "Rcpp::compileAttributes('$@')"
	@echo "==> Running oxygenize for $@"
	Rscript -e "roxygen2::roxygenize('$@', load = 'source')"
	@echo "==> making package $@"
	R CMD build --no-build-vignettes $@
