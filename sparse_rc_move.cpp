// Move-constructing CppAD::sparse_rc warns on GCC at -O2 -Wall.
// sparse_rc(sparse_rc&&) calls swap before nr_, nc_, and nnz_ are initialized,
// and std::swap reads them (bits/move.h). Including Eigen and moving the
// pattern as a member, as below, is enough for GCC 14 through 16 to report it.
//
//   g++ sparse_rc_move.cpp -o sparse_rc_move -std=gnu++20 -O2 -Wall -DNDEBUG \
//       -Werror=uninitialized \
//       -I RCppAD/inst/include \
//       -I "$(Rscript -e 'cat(system.file("include", package="RcppEigen"))')"

#include <cppad/utility/vector.hpp>
#include <cppad/utility/sparse_rc.hpp>
#include <Eigen/Sparse>
#include <iostream>
#include <utility>

struct Holder {
    Eigen::SparseMatrix<double> H;
    CppAD::sparse_rc<CppAD::vector<size_t>> pattern;
};

struct Shard {
    Holder pack;
    explicit Shard(Holder&& p) : pack(std::move(p)) {}
};

int main() {
    Holder in;
    Shard shard(std::move(in));
    std::cout << shard.pack.H.rows() << " " << shard.pack.pattern.nr() << "\n";
    return 0;
}
