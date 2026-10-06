#! /usr/bin/env bash
# Compile sparse_rc_move.cpp, which move-constructs CppAD::sparse_rc next to
# an Eigen matrix. The checked-out sparse_rc move constructor already
# initializes nr_/nc_/nnz_; this script removes that initializer so the
# warning still appears. No git command: inside the gcc image the checkout
# is owned by another user and git archive is rejected.
# https://github.com/coin-or/CppAD/issues/259
set -e -u
# -----------------------------------------------------------------------------
root="${GITHUB_WORKSPACE:-}"
if [ -z "$root" ]; then
  root=$(cd "$(dirname "$0")/../.." && pwd)
fi
src="$root/sparse_rc_move.cpp"
if [ ! -f "$src" ]; then
  echo "missing $src" >&2
  exit 1
fi
# -----------------------------------------------------------------------------
# Upstream Eigen headers (header-only).
eigen_ver=3.4.0
eigen_inc="$PWD/eigen-${eigen_ver}"
if [ ! -f "$eigen_inc/Eigen/Sparse" ]; then
  curl -fsSL -o eigen.tar.gz \
    "https://gitlab.com/libeigen/eigen/-/archive/${eigen_ver}/eigen-${eigen_ver}.tar.gz"
  tar -xzf eigen.tar.gz
fi
if [ ! -f "$eigen_inc/Eigen/Sparse" ]; then
  echo "Eigen ${eigen_ver} headers not found at ${eigen_inc}" >&2
  exit 1
fi
# -----------------------------------------------------------------------------
# Copy the vendored CppAD headers and drop the move-ctor initializer only.
# The default constructor uses the same initializer and must keep it.
cppad_inc="$PWD/cppad-include"
rm -rf "$cppad_inc"
mkdir -p "$cppad_inc"
cp -a "$root/RCppAD/inst/include/." "$cppad_inc/"
hdr="$cppad_inc/cppad/utility/sparse_rc.hpp"
awk '
  /sparse_rc\(sparse_rc&& other\)/ {
    print
    getline
    if ($0 ~ /nr_\(0\), nc_\(0\), nnz_\(0\)/) next
    print
    next
  }
  { print }
' "$hdr" > "$hdr.unfixed"
mv "$hdr.unfixed" "$hdr"
move_next=$(grep -A1 'sparse_rc(sparse_rc&& other)' "$hdr" | tail -n 1)
case "$move_next" in
  *'nr_(0), nc_(0), nnz_(0)'*)
    echo "sparse_rc move constructor still has the initializer" >&2
    exit 1
    ;;
esac
# -----------------------------------------------------------------------------
# Eigen's -Wignored-attributes notes are not this bug and do not fail the build.
cmd=(
  g++ "$src" -o sparse_rc_move
  -std=gnu++20
  -O2
  -Wall
  -DNDEBUG
  -Werror=uninitialized
  -I"$cppad_inc"
  -I"$eigen_inc"
)
printf '%q ' "${cmd[@]}"
echo
"${cmd[@]}"
./sparse_rc_move
