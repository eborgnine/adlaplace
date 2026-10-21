# RCppAD patches on top of CppAD 20260000.0

Applied by `tools/update-cppad.sh` (re-run after each upstream refresh):

1. `configure.hpp` — from `tools/configure.hpp.in` (ColPack/Eigen/ADOLC/IPOPT off).
2. `utility/error_handler.hpp` — `REprintf` + `Rf_error` (no `cerr`/`exit`).
3. `local/temp_file.hpp` — inline definition (headers-only; no libcppad_lib).
4. `local/var_op/atomic_op.hpp`, `core/ad_fun.hpp`, `core/fun_construct.hpp` —
   `std::cout` → `RCppAD::cppad_trace_stream()`.
5. `utility/sparse_rc.hpp` — move ctor initializes `nr_`/`nc_`/`nnz_` before
   `swap` (avoids GCC 14 `-Wuninitialized` on Windows).
6. `core/fun_construct.hpp` — `ADFun(ADFun&&)` delegates to `ADFun()` before
   `swap` (avoids an UndefinedBehaviorSanitizer invalid-`bool` load).

Helper: `inst/include/RCppAD/cppad_trace_stream.hpp`.
