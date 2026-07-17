#ifndef RCPPAD_CPPAD_TRACE_STREAM_HPP
#define RCPPAD_CPPAD_TRACE_STREAM_HPP

/*
  Discard ostream used where vendored CppAD headers would otherwise
  reference std::cout (CRAN "compiled code" NOTE). Trace paths in CppAD
  are off by default; this sink exists so those branches do not pull
  std::cout into dependent shared libraries.
*/

#include <ostream>
#include <streambuf>

namespace RCppAD {

struct cppad_null_streambuf : std::streambuf {
  int overflow(int c) override { return c; }
};

inline std::ostream &cppad_trace_stream() {
  static cppad_null_streambuf buf;
  static std::ostream os(&buf);
  return os;
}

} // namespace RCppAD

#endif
