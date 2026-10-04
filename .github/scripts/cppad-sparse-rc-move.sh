#! /usr/bin/env bash
# Compile sparse_rc_move.cpp, which move-constructs CppAD::sparse_rc next to
# an Eigen matrix. Headers are CppAD from b804867, before the move constructor
# initialized nr_/nc_/nnz_. The copy in this repo already has that fix.
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
# Pre-fix CppAD include tree (sparse_rc move ctor does not initialize scalars).
cppad_commit=b80486726a580fbd2820a3dbed73c183bbc3c526
cppad_inc="$PWD/cppad-include"
rm -rf "$cppad_inc"
mkdir -p "$cppad_inc"
git -C "$root" archive "$cppad_commit" RCppAD/inst/include | tar -x -C "$cppad_inc"
cppad_inc="$cppad_inc/RCppAD/inst/include"
move_next=$(grep -A1 'sparse_rc(sparse_rc&& other)' "$cppad_inc/cppad/utility/sparse_rc.hpp" | tail -n 1)
case "$move_next" in
  *'nr_(0), nc_(0), nnz_(0)'*)
    echo "pre-fix sparse_rc still has the move-ctor initializer" >&2
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
