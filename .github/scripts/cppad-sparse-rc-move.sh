#! /usr/bin/env bash
# Minimal stand-in for the adlaplace -Wuninitialized report on Rtools.
# The move constructor reads nr_/nc_/nnz_ before they are set, then a holder
# that also contains an Eigen::SparseMatrix is moved, as AdTape is.
# The Rtools diagnostic shows up in a TU that includes Eigen; the three-scalar
# version without Eigen did not warn.
# https://github.com/coin-or/CppAD/issues/259
set -e -u
# -----------------------------------------------------------------------------
# Upstream Eigen headers (header-only). Not the copy shipped by RcppEigen.
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
#
# temp.cpp
cat << EOF > temp.cpp
#include <Eigen/Sparse>
#include <cstddef>
#include <iostream>
#include <utility>
//
// Move ctor swaps nr_/nc_/nnz_ before they are initialized.
class sparse_rc {
public:
    std::size_t nr_;
    std::size_t nc_;
    std::size_t nnz_;
    sparse_rc(void)
    : nr_(0), nc_(0), nnz_(0)
    { }
    sparse_rc(sparse_rc&& other)
    {  swap(other); }
    void swap(sparse_rc& other)
    {  std::swap(nr_, other.nr_);
       std::swap(nc_, other.nc_);
       std::swap(nnz_, other.nnz_);
    }
};
//
// Eigen matrix in the same object. Moving it is what makes Rtools warn.
struct Holder {
    Eigen::SparseMatrix<double> H;
    sparse_rc pattern;
};
struct Shard {
    Holder pack;
    explicit Shard(Holder&& p)
    : pack(std::move(p))
    { }
};
//
int main(void)
{  Holder in;
   Shard shard(std::move(in));
   std::cout << shard.pack.H.rows() << " " << shard.pack.pattern.nr_ << "\n";
   return 0;
}
EOF
#
# Rtools adlaplace line, plus -Werror=uninitialized.
# Eigen's -Wignored-attributes notes are not this bug and do not fail the build.
cmd=(
  g++ temp.cpp -o temp
  -std=gnu++20
  -O2
  -Wall
  -DNDEBUG
  -Werror=uninitialized
  -I"$eigen_inc"
)
printf '%q ' "${cmd[@]}"
echo
"${cmd[@]}"
./temp
