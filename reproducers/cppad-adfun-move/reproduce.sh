#!/usr/bin/env bash
# Reproduce the CppAD 20260000.0 ADFun move-constructor UBSAN error.
# No R. Downloads the upstream tag, compiles adfun_move.cpp with
# -fsanitize=undefined, and succeeds only when the invalid-bool load appears.
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
ver="${CPPAD_VERSION:-20260000.0}"
cxx="${CXX:-c++}"
cache="${TMPDIR:-/tmp}/cppad-${ver}"

if [ ! -f "$cache/include/cppad/core/fun_construct.hpp" ]; then
  mkdir -p "$cache"
  curl -fsSL "https://github.com/coin-or/CppAD/archive/refs/tags/${ver}.tar.gz" \
    | tar -xz -C "$cache" --strip-components=1
fi

work="$(mktemp -d "${TMPDIR:-/tmp}/cppad-adfun-ubsan.XXXXXX")"
trap 'rm -rf "$work"' EXIT
cp -a "$cache/include/." "$work/include"
# The tag ships configure.hpp.in. This generated header does not change the move constructor.
cp "$root/configure.hpp" "$work/include/cppad/configure.hpp"

"$cxx" -std=c++17 -O1 -g \
  -fsanitize=undefined -fno-sanitize-recover=undefined \
  -I"$work/include" \
  "$root/adfun_move.cpp" -o "$work/adfun_move"

set +e
"$work/adfun_move" >"$work/log" 2>&1
status=$?
set -e
cat "$work/log"

if [ "$status" -ne 0 ] && grep -q "not a valid value for type" "$work/log"; then
  echo "reproduced: CppAD ${ver} ADFun move reads an invalid bool (exit ${status})"
  exit 0
fi

echo "did not reproduce the CppAD ADFun UBSAN error (exit ${status})" >&2
exit 1
