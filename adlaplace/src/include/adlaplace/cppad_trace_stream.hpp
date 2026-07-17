#ifndef ADLAPLACE_CPPAD_TRACE_STREAM_HPP
#define ADLAPLACE_CPPAD_TRACE_STREAM_HPP

/*
  Discard ostream used by configure-generated CppAD header shadows that
  replace std::cout (CRAN "compiled code" NOTE). Trace paths in CppAD are
  off by default; this sink exists only so those branches do not reference
  std::cout in the shared library.
*/

#include <ostream>
#include <streambuf>

namespace adlaplace {

struct cppad_null_streambuf : std::streambuf {
  int overflow(int c) override { return c; }
};

inline std::ostream &cppad_trace_stream() {
  static cppad_null_streambuf buf;
  static std::ostream os(&buf);
  return os;
}

} // namespace adlaplace

#endif
