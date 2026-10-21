// Minimal CppAD ADFun move-constructor reproducer.
//
// CppAD::ADFun<Base>::ADFun(ADFun&&) calls swap() before its scalar members
// exist. std::swap then loads those bools. Filling the destination storage
// with 0xBE makes the load deterministic under -fsanitize=undefined:
// value 190 is not a valid bool.
//
// Build and run with reproducers/cppad-adfun-move/run.sh.

#include <cppad/cppad.hpp>

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <new>
#include <utility>

int main() {
  using Fun = CppAD::ADFun<double>;

  Fun source;

  alignas(Fun) unsigned char storage[sizeof(Fun)];
  std::memset(storage, 0xBE, sizeof(storage));

  Fun* destination = new (storage) Fun(std::move(source));
  const std::size_t n_var = destination->size_var();
  destination->~Fun();

  std::printf("moved ADFun size_var=%zu\n", n_var);
  return n_var == 0 ? 0 : 1;
}
