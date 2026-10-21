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
