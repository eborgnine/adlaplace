# RCppAD 20260000.0-3

* Document that the `ADFun(ADFun&&)` initializer matches upstream
  coin-or/CppAD@dbf793e (2026-11-09). The constructor form is unchanged from
  1.20260000.0-2. A full refresh to CppAD 20261109 is not included: that tree
  moves `thread_alloc` into `cppad_lib`, which is incompatible with this
  headers-only package.

# RCppAD 20260000.0-2

* `ADFun(ADFun&&)` now delegates to the default constructor before `swap`.
  The previous constructor read uninitialized `bool` members, which
  UndefinedBehaviorSanitizer reports while `adlaplace` move-constructs tapes.

# RCppAD 20260000.0-1

* Initial CRAN-oriented release: vendored CppAD 20260000.0 headers with
  ColPack/Eigen/ADOLC/IPOPT disabled in `configure.hpp`, R-safe
  `ErrorHandler` (`REprintf` / `Rf_error`), and `std::cout` redirects to a
  discard stream for CRAN compiled-code checks.
* Inline `CppAD::local::temp_file()` so LinkingTo packages do not need
  `libcppad_lib` (headers-only packaging).
