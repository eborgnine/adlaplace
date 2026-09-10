#ifndef ADLAPLACE_EIGEN_DIAGNOSTICS_HPP
#define ADLAPLACE_EIGEN_DIAGNOSTICS_HPP

// Wrap Eigen / RcppEigen includes with these macros so Rtools GCC 14
// false-positive -Wuninitialized (std::move of packet types in bits/move.h)
// does not become a CRAN "significant warning". Do not put -Wno-* on
// PKG_CXXFLAGS (winbuilder treats that as a non-portable flag WARNING).
//
// Clang understands GCC diagnostic but not -Wmaybe-uninitialized; enabling
// that pragma on Clang can itself become a significant warning.

#ifdef __GNUC__
# define ADLAPLACE_EIGEN_DIAG_PUSH \
  _Pragma("GCC diagnostic push") \
  _Pragma("GCC diagnostic ignored \"-Wuninitialized\"")
# ifndef __clang__
#  define ADLAPLACE_EIGEN_DIAG_PUSH_MAYBE \
     ADLAPLACE_EIGEN_DIAG_PUSH \
     _Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# else
#  define ADLAPLACE_EIGEN_DIAG_PUSH_MAYBE ADLAPLACE_EIGEN_DIAG_PUSH
# endif
# define ADLAPLACE_EIGEN_DIAG_POP _Pragma("GCC diagnostic pop")
#else
# define ADLAPLACE_EIGEN_DIAG_PUSH_MAYBE
# define ADLAPLACE_EIGEN_DIAG_POP
#endif

#endif
