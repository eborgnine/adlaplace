# Build and install monorepo packages in dependency order.
#
#   make            # document, build, and install everything
#   make adlaplace  # one package (deps are installed first)
#   make check      # R CMD check each package's .tar.gz (run make first)
#   make crancheck  # same with --as-cran
#   make clean      # cleanup scripts + src/*.o,src/*.so
#
# Order: RCppAD → adlaplace → {Example, Hgp, Fem, admvn} → hpolcc
# (hpolcc needs adlaplaceHgp; backends need adlaplace headers.)

INSTALL_ORDER := RCppAD adlaplace adlaplaceExample adlaplaceHgp adlaplaceFem hpolcc

.DEFAULT_GOAL := all

.PHONY: all install check crancheck clean $(INSTALL_ORDER) dirichlet_multinom.pdf

# Avoid concurrent R CMD INSTALL locking the same library when using make -j.
.NOTPARALLEL: $(INSTALL_ORDER)

all: install

install: $(INSTALL_ORDER)

# --- dependency edges -------------------------------------------------------

adlaplace: RCppAD
adlaplaceExample: adlaplace
adlaplaceHgp: adlaplace
adlaplaceFem: adlaplace
admvn: RCppAD
hpolcc: adlaplace adlaplaceHgp

# --- per-package: cleanup → attributes → roxygen → build → install ----------

$(INSTALL_ORDER):
	@echo ""
	@echo "==> $@"
	@if [ -f "$@/cleanup" ]; then \
		echo "==> cleanup"; \
		(cd "$@" && sh ./cleanup); \
	fi
	@if [ -d "$@/src" ] && ls "$@"/src/*.cpp >/dev/null 2>&1; then \
		echo "==> Rcpp::compileAttributes"; \
		Rscript -e "Rcpp::compileAttributes('$@')"; \
	fi
	@echo "==> roxygen2::roxygenize"
	@Rscript -e "roxygen2::roxygenize('$@')"
	@echo "==> R CMD build"
	@R CMD build --no-build-vignettes "$@"
	@PKG="$$(sed -n 's/^Package:[[:space:]]*//p' "$@/DESCRIPTION" | head -n 1)"; \
	VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' "$@/DESCRIPTION" | head -n 1)"; \
	TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
	test -f "$$TARBALL" || { echo "Expected tarball $$TARBALL not found"; exit 1; }; \
	echo "==> R CMD INSTALL $$TARBALL"; \
	R CMD INSTALL "$$TARBALL"

# --- check all built tarballs -----------------------------------------------

define require_tarballs
	missing=0; \
	for pkg in $(INSTALL_ORDER); do \
		PKG="$$(sed -n 's/^Package:[[:space:]]*//p' "$$pkg/DESCRIPTION" | head -n 1)"; \
		VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' "$$pkg/DESCRIPTION" | head -n 1)"; \
		TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
		if [ ! -f "$$TARBALL" ]; then \
			echo "Missing $$TARBALL (run 'make' first)"; \
			missing=1; \
		fi; \
	done; \
	if [ "$$missing" -ne 0 ]; then exit 1; fi
endef

define check_tarballs
	for pkg in $(INSTALL_ORDER); do \
		PKG="$$(sed -n 's/^Package:[[:space:]]*//p' "$$pkg/DESCRIPTION" | head -n 1)"; \
		VERSION="$$(sed -n 's/^Version:[[:space:]]*//p' "$$pkg/DESCRIPTION" | head -n 1)"; \
		TARBALL="$${PKG}_$${VERSION}.tar.gz"; \
		echo ""; \
		echo "==> R CMD check $(1) $$TARBALL"; \
		R CMD check $(1) "$$TARBALL" || exit 1; \
	done
endef

check:
	@$(require_tarballs)
	@$(call check_tarballs,)

crancheck:
	@$(require_tarballs)
	@$(call check_tarballs,--as-cran)

# --- housekeeping -----------------------------------------------------------

clean:
	@echo "==> Cleaning package build artifacts (cleanup scripts + src objects)"
	@for pkg in $(INSTALL_ORDER); do \
		if [ -f "$$pkg/cleanup" ]; then \
			echo "==> $$pkg/cleanup"; \
			(cd "$$pkg" && sh ./cleanup); \
		fi; \
	done
	@find . -type f \( -name '*.o' -o -name '*.so' \) -path '*/src/*' -delete
	@rm -f $(foreach p,$(INSTALL_ORDER),$(p)_*.tar.gz)

dirichlet_multinom.pdf: adlaplace/vignettes/dirichlet_multinom.Rmd
	pandoc adlaplace/vignettes/dirichlet_multinom.Rmd -o dirichlet_multinom.pdf
