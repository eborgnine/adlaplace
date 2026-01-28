
# Make the default target run all packages
.DEFAULT_GOAL := all

# discover package dirs (subdirectories containing DESCRIPTION)
PKGS := $(patsubst %/DESCRIPTION,%,$(wildcard */DESCRIPTION))

.PHONY: all $(PKGS)

all: $(PKGS)

$(PKGS):
	@echo "==> Running compileAttributes for $@"
	Rscript -e "Rcpp::compileAttributes('$@')"
	@echo "==> Running roxygenize for $@"
	Rscript -e "roxygen2::roxygenize('$@', load = 'source')"
	@echo "==> making package $@"
	R CMD build --no-build-vignettes $@