#!/bin/sh
# Refresh vendored CppAD headers into inst/include/cppad/.
# Maintainer tool only (listed in .Rbuildignore). Not run at R CMD INSTALL.
#
# Usage:
#   ./tools/update-cppad.sh 20260000.0
#   ./tools/update-cppad.sh /path/to/include   # directory containing cppad/
#   ./tools/update-cppad.sh                   # try Homebrew Cellar

set -eu

ROOT="$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)"
TOOLS="$ROOT/tools"
DEST="$ROOT/inst/include/cppad"
VERSION_ARG="${1:-}"

resolve_src () {
  if [ -n "$VERSION_ARG" ] && [ -d "$VERSION_ARG/cppad" ]; then
    printf '%s' "$VERSION_ARG"
    return 0
  fi
  if [ -n "$VERSION_ARG" ] && [ -d "$VERSION_ARG/include/cppad" ]; then
    printf '%s' "$VERSION_ARG/include"
    return 0
  fi
  if [ -n "$VERSION_ARG" ]; then
    # Version tag: try Homebrew, then download from GitHub
    ver="$VERSION_ARG"
    for base in \
      "/opt/homebrew/Cellar/cppad/$ver/include" \
      "/usr/local/Cellar/cppad/$ver/include"
    do
      if [ -d "$base/cppad" ]; then
        printf '%s' "$base"
        return 0
      fi
    done
    echo "update-cppad: downloading CppAD $ver from GitHub..." >&2
    tmp="$(mktemp -d "${TMPDIR:-/tmp}/rcppad-cppad.XXXXXX")"
    url="https://github.com/coin-or/CppAD/archive/refs/tags/$ver.tar.gz"
    if ! curl -fsSL "$url" | tar -xz -C "$tmp"; then
      echo "update-cppad: ERROR: failed to fetch $url" >&2
      rm -rf "$tmp"
      return 1
    fi
    # Archive top dir is usually CppAD-<ver>
    src_dir="$(find "$tmp" -maxdepth 2 -type d -name include | head -n 1)"
    if [ -z "$src_dir" ] || [ ! -d "$src_dir/cppad" ]; then
      echo "update-cppad: ERROR: no include/cppad in downloaded archive" >&2
      rm -rf "$tmp"
      return 1
    fi
    printf '%s' "$src_dir"
    echo "$tmp" > "$tmp/.rcppad_tmpdir"
    return 0
  fi
  # Default: newest Homebrew cppad
  for cellar in /opt/homebrew/Cellar/cppad /usr/local/Cellar/cppad; do
    if [ -d "$cellar" ]; then
      latest="$(ls -1 "$cellar" | sort -V | tail -n 1)"
      if [ -n "$latest" ] && [ -d "$cellar/$latest/include/cppad" ]; then
        printf '%s' "$cellar/$latest/include"
        return 0
      fi
    fi
  done
  echo "update-cppad: ERROR: pass a version, include path, or install Homebrew cppad" >&2
  return 1
}

SRC_INCLUDE="$(resolve_src)"
CLEANUP_TMP=""
if [ -f "${SRC_INCLUDE%/*}/.rcppad_tmpdir" ] 2>/dev/null; then
  CLEANUP_TMP="$(cat "${SRC_INCLUDE%/*}/.rcppad_tmpdir" 2>/dev/null || true)"
fi
# When we downloaded, tmp marker is next to archive root
parent_of_include="$(dirname "$SRC_INCLUDE")"
if [ -f "$parent_of_include/.rcppad_tmpdir" ]; then
  CLEANUP_TMP="$(cat "$parent_of_include/.rcppad_tmpdir")"
fi
# Downloaded layout: tmp/CppAD-ver/include — marker on tmp
grand="$(dirname "$parent_of_include")"
if [ -f "$grand/.rcppad_tmpdir" ]; then
  CLEANUP_TMP="$(cat "$grand/.rcppad_tmpdir")"
fi

if [ ! -f "$SRC_INCLUDE/cppad/cppad.hpp" ]; then
  echo "update-cppad: ERROR: missing $SRC_INCLUDE/cppad/cppad.hpp" >&2
  exit 1
fi

# Infer version from PACKAGE_STRING in upstream configure if present
VERSION="$VERSION_ARG"
if [ -z "$VERSION" ] || [ -d "$VERSION_ARG" ] || [ -d "${VERSION_ARG:-}/cppad" ]; then
  if [ -f "$SRC_INCLUDE/cppad/configure.hpp" ]; then
    VERSION="$(sed -n 's/.*CPPAD_PACKAGE_STRING "cppad-\([^"]*\)".*/\1/p' \
      "$SRC_INCLUDE/cppad/configure.hpp" | head -n 1)"
  fi
fi
if [ -z "$VERSION" ]; then
  VERSION="unknown"
fi

echo "update-cppad: source=$SRC_INCLUDE/cppad version=$VERSION" >&2

rm -rf "$DEST"
mkdir -p "$(dirname "$DEST")"
cp -a "$SRC_INCLUDE/cppad" "$DEST"

# Frozen configure.hpp
sed "s/@CPPAD_VERSION@/$VERSION/g" "$TOOLS/configure.hpp.in" > "$DEST/configure.hpp"

# R-safe ErrorHandler (replace upstream utility/error_handler.hpp)
cp "$TOOLS/patches/error_handler.hpp" "$DEST/utility/error_handler.hpp"

# Inline temp_file() so LinkingTo packages need no libcppad_lib
cp "$TOOLS/patches/temp_file.hpp" "$DEST/local/temp_file.hpp"

# Redirect std::cout in headers that appear in CRAN compiled-code checks
redirect_cout () {
  rel="$1"
  f="$DEST/$rel"
  if [ ! -f "$f" ]; then
    echo "update-cppad: WARNING: missing $rel (skip cout redirect)" >&2
    return 0
  fi
  tmp="$f.rcppad.tmp"
  {
    printf '%s\n' \
      "/* RCppAD: std::cout -> RCppAD::cppad_trace_stream() (CRAN compiled code). */" \
      "#include <RCppAD/cppad_trace_stream.hpp>"
    sed -e 's/std::cout/::RCppAD::cppad_trace_stream()/g' "$f"
  } > "$tmp"
  mv "$tmp" "$f"
}

redirect_cout "local/var_op/atomic_op.hpp"
redirect_cout "core/ad_fun.hpp"
redirect_cout "core/fun_construct.hpp"

# Document patches
cat > "$DEST/README.RCppAD.md" <<EOF
# RCppAD patches on top of CppAD $VERSION

Applied by \`tools/update-cppad.sh\` (re-run after each upstream refresh):

1. \`configure.hpp\` — from \`tools/configure.hpp.in\` (ColPack/Eigen/ADOLC/IPOPT off).
2. \`utility/error_handler.hpp\` — \`REprintf\` + \`Rf_error\` (no \`cerr\`/\`exit\`).
3. \`local/temp_file.hpp\` — inline definition (headers-only; no libcppad_lib).
4. \`local/var_op/atomic_op.hpp\`, \`core/ad_fun.hpp\`, \`core/fun_construct.hpp\` —
   \`std::cout\` → \`RCppAD::cppad_trace_stream()\`.

Helper: \`inst/include/RCppAD/cppad_trace_stream.hpp\`.
EOF

# Bump DESCRIPTION Version / Date
desc="$ROOT/DESCRIPTION"
pkg_ver="${VERSION}-1"
if command -v sed >/dev/null 2>&1; then
  tmpdesc="$desc.tmp"
  sed -e "s/^Version:.*/Version: $pkg_ver/" \
      -e "s/^Date:.*/Date: $(date +%Y-%m-%d)/" \
      "$desc" > "$tmpdesc" || cp "$desc" "$tmpdesc"
  # Date field may be absent
  if ! grep -q '^Date:' "$tmpdesc"; then
    sed -e "s/^Version:.*/Version: $pkg_ver/" "$desc" > "$tmpdesc"
  else
    sed -e "s/^Version:.*/Version: $pkg_ver/" \
        -e "s/^Date:.*/Date: $(date +%Y-%m-%d)/" \
        "$desc" > "$tmpdesc"
  fi
  mv "$tmpdesc" "$desc"
fi

nfiles="$(find "$DEST" -type f | wc -l | tr -d ' ')"
echo "update-cppad: wrote $nfiles files under inst/include/cppad/" >&2
echo "update-cppad: DESCRIPTION Version -> $pkg_ver" >&2

if [ -n "$CLEANUP_TMP" ] && [ -d "$CLEANUP_TMP" ]; then
  rm -rf "$CLEANUP_TMP"
fi

exit 0
