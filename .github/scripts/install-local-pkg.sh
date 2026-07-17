#!/usr/bin/env bash
# Build a source tarball and install from it (not from the raw package directory).
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: $0 <package-dir>" >&2
  exit 1
fi

pkg=$1
if [ ! -d "$pkg" ]; then
  echo "error: package directory not found: $pkg" >&2
  exit 1
fi

# Windows checkout may omit +x even when git records 100755.
chmod +x "$pkg/configure" "$pkg/cleanup" "$pkg/configure.win" "$pkg/cleanup.win" 2>/dev/null || true

R CMD build "$pkg" --no-manual --no-build-vignettes
tarball=$(ls -t "${pkg}"_*.tar.gz | head -n 1)
if [ -z "$tarball" ] || [ ! -f "$tarball" ]; then
  echo "error: built tarball not found for $pkg" >&2
  exit 1
fi

R CMD INSTALL "$tarball"
rm -f "$tarball"
