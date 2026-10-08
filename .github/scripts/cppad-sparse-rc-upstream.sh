#! /usr/bin/env bash
# Compile sparse_rc_move.cpp against upstream CppAD commit 302e241
# (coin-or/CppAD#259). That commit initializes nr_/nc_/nnz_ in the
# sparse_rc move constructor. This script does not use the vendored
# RCppAD headers and does not strip the initializer.
# https://github.com/coin-or/CppAD/commit/302e241ae7b6b146c676590981b5919af71e1ccf
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
# Upstream CppAD at the issue 259 fix. Not a release tag.
cppad_sha=302e241ae7b6b146c676590981b5919af71e1ccf
curl -fsSL -o cppad.tar.gz \
  "https://github.com/coin-or/CppAD/archive/${cppad_sha}.tar.gz"
tar -xzf cppad.tar.gz
cppad_hpp="$(find "$PWD" -type f -path "*/include/cppad/cppad.hpp" -print | head -n 1)"
if [ -z "$cppad_hpp" ]; then
  echo "CppAD ${cppad_sha} headers not found" >&2
  exit 1
fi
# -I needs the include/ directory, not include/cppad/.
cppad_inc="$(dirname "$(dirname "$cppad_hpp")")"
# A git archive ships configure.hpp.in. CMake writes configure.hpp, and
# vector.hpp includes it. Fill the upstream template with header-only
# defaults. This does not change sparse_rc.hpp.
cfg_in="$cppad_inc/cppad/configure.hpp.in"
cfg_out="$cppad_inc/cppad/configure.hpp"
if [ ! -f "$cfg_in" ]; then
  echo "missing $cfg_in" >&2
  exit 1
fi
cppad_ver="$(sed -n 's/.*SET(cppad_version "\([^"]*\)").*/\1/p' \
  "$cppad_inc/../CMakeLists.txt" | head -n 1)"
if [ -z "$cppad_ver" ]; then
  echo "could not read cppad_version from CMakeLists.txt" >&2
  exit 1
fi
sed \
  -e 's/@cppad_lib_static_01@/0/g' \
  -e 's/@cppad_link_flags_has_m32@/0/g' \
  -e 's/@compiler_has_conversion_warn@/0/g' \
  -e 's/@cppad_debug_and_release_01@/1/g' \
  -e 's/@use_cplusplus_2017_ok@/1/g' \
  -e "s/@cppad_version@/${cppad_ver}/g" \
  -e 's/@cppad_has_adolc@/0/g' \
  -e 's/@cppad_has_colpack@/0/g' \
  -e 's/@cppad_has_eigen@/0/g' \
  -e 's/@cppad_has_ipopt@/0/g' \
  -e 's/@cppad_deprecated_01@/0/g' \
  -e 's/@cppad_boostvector@/0/g' \
  -e 's/@cppad_cppadvector@/1/g' \
  -e 's/@cppad_stdvector@/0/g' \
  -e 's/@cppad_eigenvector@/0/g' \
  -e 's/@cppad_has_gettimeofday@/0/g' \
  -e 's/@cppad_tape_addr_type@/unsigned int/g' \
  -e 's/@cppad_is_same_tape_addr_type_size_t@/0/g' \
  -e 's/@cppad_tape_id_type@/unsigned int/g' \
  -e 's/@cppad_max_num_threads@/48/g' \
  -e 's/@cppad_has_mkstemp@/0/g' \
  -e 's/@cppad_has_tmpnam_s@/0/g' \
  -e 's/@cppad_c_compiler_cmd@/cc/g' \
  -e 's/@cppad_c_compiler_gnu_flags@/0/g' \
  -e 's/@cppad_c_compiler_msvc_flags@/0/g' \
  -e 's/@cppad_is_same_unsigned_int_size_t@/0/g' \
  -e 's/@cppad_padding_block_t@//g' \
  "$cfg_in" > "$cfg_out"
if grep -n '^# *define' "$cfg_out" | grep -q '@'; then
  echo "configure.hpp still has unsubstituted placeholders" >&2
  exit 1
fi
hdr="$cppad_inc/cppad/utility/sparse_rc.hpp"
if [ ! -f "$hdr" ]; then
  echo "missing $hdr" >&2
  exit 1
fi
if ! grep -q 'sparse_rc(sparse_rc&& other)' "$hdr"; then
  echo "sparse_rc move constructor not found in $hdr" >&2
  exit 1
fi
move_next=$(grep -A1 'sparse_rc(sparse_rc&& other)' "$hdr" | tail -n 1)
case "$move_next" in
  *'nr_(0), nc_(0), nnz_(0)'*)
    ;;
  *)
    echo "upstream sparse_rc move constructor is missing nr_(0), nc_(0), nnz_(0)" >&2
    echo "next line was: $move_next" >&2
    exit 1
    ;;
esac
# -----------------------------------------------------------------------------
# thread_info and capacity_info live in cppad_lib as of this commit.
# CppAD::vector, used as sparse_rc's index vector, calls them.
# CPPAD_LIB_EXPORTS keeps those symbols from being dllimport on Windows.
cppad_root="$(dirname "$cppad_inc")"
thread_alloc_src="$cppad_root/cppad_lib/static/thread_alloc.cpp"
if [ ! -f "$thread_alloc_src" ]; then
  echo "missing $thread_alloc_src" >&2
  exit 1
fi
# Eigen's -Wignored-attributes notes are not this bug and do not fail the build.
cmd=(
  g++ "$src" "$thread_alloc_src" -o sparse_rc_move
  -std=gnu++20
  -O2
  -Wall
  -DNDEBUG
  -Werror=uninitialized
  -DCPPAD_LIB_EXPORTS
  -I"$cppad_inc"
  -I"$eigen_inc"
)
printf '%q ' "${cmd[@]}"
echo
"${cmd[@]}"
out="$(./sparse_rc_move | tr -d '\r')"
printf '%s\n' "$out"
case "$out" in
  '0 0')
    ;;
  *)
    echo "expected output '0 0', got '${out}'" >&2
    exit 1
    ;;
esac
