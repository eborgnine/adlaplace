#! /usr/bin/env bash
# Minimal stand-in for the adlaplace -Wuninitialized report.
# Unlike Brad's swap.sh, the uninitialized scalars are read by the move
# constructor itself, then that constructor runs as part of moving an
# outer struct (the AdTape chain in the Rtools log).
# https://github.com/coin-or/CppAD/issues/259
set -e -u
# -----------------------------------------------------------------------------
echo_eval() {
   echo $*
   eval $*
}
# -----------------------------------------------------------------------------
#
# temp.cpp
cat << EOF > temp.cpp
#include <cstddef>
#include <iostream>
#include <utility>
//
// sparse_rc: move ctor has no initializer and swap reads nr_/nc_/nnz_.
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
// AdTape: implicit move constructor move-constructs unused_pattern.
struct AdTape {
    sparse_rc unused_pattern;
};
//
// Passing and returning by value both call AdTape's move constructor.
AdTape relocate(AdTape src)
{  return src; }
//
int main(void)
{  AdTape in;
   AdTape out = relocate(std::move(in));
   std::cout << out.unused_pattern.nr_ << "\n";
   return 0;
}
EOF
#
# Same flags as Brad's swap.sh.
echo_eval g++ temp.cpp -o temp \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -Wshadow \
    -Wconversion \
    -Wlogical-op \
    -Wduplicated-cond \
    -Wduplicated-branches \
    -Wunused \
    -Wold-style-cast \
    -Woverloaded-virtual \
    -Wnull-dereference \
    -Wformat=2 -Werror
#
echo_eval ./temp
