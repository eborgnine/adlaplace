# CppAD `ADFun` move constructor reads an uninitialized `bool`

CppAD 20260000.0 (`ADFun(ADFun&&)`) calls `swap()` before the destination object exists. `std::swap` then loads `exceed_collision_limit_`, `has_been_optimized_`, and `check_for_nan_`. UndefinedBehaviorSanitizer reports that load when the destination storage is filled with `0xBE` (decimal 190), which is not a valid `bool`.

This does not depend on R, Eigen, or OpenMP. It is the same sanitizer finding CRAN reports while `adlaplace` move-constructs an `AdTape` into a `std::vector`:

- [macOS arm64](https://www.stats.ox.ac.uk/pub/bdr/M1-SAN/adlaplace/tests/testthat.Rout)
- [GCC UBSAN](https://www.stats.ox.ac.uk/pub/bdr/memtests/gcc-UBSAN/adlaplace/tests/testthat.Rout)

## Fault

`include/cppad/core/fun_construct.hpp`:

```cpp
template <class Base, class RecBase>
ADFun<Base,RecBase>::ADFun(ADFun&& f)
{  swap(f); }
```

`swap()` exchanges the three `bool` members with `std::swap`. The move constructor never initializes them first. The comment above this constructor says those default values do not matter to the destructor; they do matter to `swap()`.

## Reproduce

The tag archive ships `configure.hpp.in` rather than `configure.hpp`. `run.sh --fetch` downloads CppAD 20260000.0 and overlays `configure.hpp` from this directory. That generated header does not change `fun_construct.hpp`.

```bash
./run.sh --fetch 20260000.0 --cxx clang++
./run.sh --fetch 20260000.0 --cxx g++
```

Both use:

```text
-std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined
```

Apple Clang 21 exits 134. The full log is [upstream-ubsan.txt](upstream-ubsan.txt):

```text
runtime error: load of value 190, which is not a valid value for type 'bool'
SUMMARY: UndefinedBehaviorSanitizer: undefined-behavior .../__utility/swap.h:43:11
```

GCC's libstdc++ reports the same load from `bits/move.h`, as in the Oxford GCC log linked above.

## Fix

Delegate to the default constructor, which initializes the scalars, and then swap:

```cpp
template <class Base, class RecBase>
ADFun<Base,RecBase>::ADFun(ADFun&& f)
: ADFun<Base, RecBase>()
{  swap(f); }
```

```bash
./run.sh --fetch 20260000.0 --patch-move --cxx clang++
```

That command exits 0 and prints `moved ADFun size_var=0`. The log is [patched-ubsan.txt](patched-ubsan.txt).

`--vendored` overlays `RCppAD`'s patched `fun_construct.hpp` instead of editing the tag with `awk`. Its log is [vendored-ubsan.txt](vendored-ubsan.txt).
