#include "pmvn_atomic.hpp"

namespace admvn {
namespace detail {

namespace {
thread_local std::vector<MvnTape*> pmvn_atomic_tapes;
}  // namespace

MvnTape* pmvn_atomic_tape_from_id(size_t call_id) {
  if (call_id >= pmvn_atomic_tapes.size()) {
    return nullptr;
  }
  return pmvn_atomic_tapes[call_id];
}

void set_pmvn_atomic_tape(size_t call_id, MvnTape* tape) {
  if (call_id >= pmvn_atomic_tapes.size()) {
    pmvn_atomic_tapes.resize(call_id + 1, nullptr);
  }
  pmvn_atomic_tapes[call_id] = tape;
}

void set_pmvn_atomic_tapes(MvnTape* p1, MvnTape* p2) {
  set_pmvn_atomic_tape(0, p1);
  set_pmvn_atomic_tape(1, p2);
}

}  // namespace detail
}  // namespace admvn
