# Force: do not use ColPack paths
CPPAD_CPPFLAGS="-DCPPAD_HAS_COLPACK=0"
CPPAD_INCLUDE_CPPFLAGS=""
CPPAD_LIBS=""

OPENMP_CPPFLAGS=""
OPENMP_CXXFLAGS=""
OPENMP_LIBS=""
BREW_CPPFLAGS=""
BREW_PREFIX=""

UNAME_S="$(uname -s 2>/dev/null || echo unknown)"
R_HOME_DIR="$(R RHOME 2>/dev/null || true)"

CXX_CMD="$(R CMD config CXX 2>/dev/null || echo c++)"
CXXFLAGS_CMD="$(R CMD config CXXFLAGS 2>/dev/null || true)"
CPPFLAGS_CMD="$(R CMD config CPPFLAGS 2>/dev/null || true)"

# Eigen: packages LinkingTo RcppEigen get -I from R; do not probe here.

# Probe scratch dirs go under TMPDIR (mktemp). Relative config.$$ in the
# package tree is fragile during R CMD INSTALL / staged installs.
make_conftest_dir () {
  mktemp -d "${TMPDIR:-/tmp}/adlaplace-conftest.XXXXXX"
}

check_cxxflag () {
  flag="$1"
  tmpdir="$(make_conftest_dir)" || return 1
  cat > "$tmpdir/conftest.cpp" <<'EOF'
int main() { return 0; }
EOF
  if eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD $flag -c $tmpdir/conftest.cpp -o $tmpdir/conftest.o" \
       >/dev/null 2>&1; then
    rm -rf "$tmpdir"
    return 0
  fi
  echo "configure: C++ flag probe failed for '$flag':" >&2
  eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD $flag -c $tmpdir/conftest.cpp -o $tmpdir/conftest.o" >&2 || true
  rm -rf "$tmpdir"
  return 1
}

check_openmp_darwin () {
  tmpdir="$(make_conftest_dir)" || return 1
  cat > "$tmpdir/conftest.cpp" <<'EOF'
#include <omp.h>
int main() {
  omp_set_num_threads(1);
  return omp_get_max_threads() >= 1 ? 0 : 1;
}
EOF
  # Prefer dynamic libomp: static libomp.a + dyn.load (pkgload/load_all) can
  # segfault in OpenMP TLS on macOS under R's -undefined dynamic_lookup.
  if eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD -I$BREW_PREFIX/opt/libomp/include -Xpreprocessor -fopenmp \
      $tmpdir/conftest.cpp -L$BREW_PREFIX/opt/libomp/lib -lomp \
      -Wl,-rpath,$BREW_PREFIX/opt/libomp/lib -o $tmpdir/conftest" \
      >/dev/null 2>&1; then
    rm -rf "$tmpdir"
    return 0
  fi
  echo "configure: macOS OpenMP link probe failed:" >&2
  eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD -I$BREW_PREFIX/opt/libomp/include -Xpreprocessor -fopenmp \
      $tmpdir/conftest.cpp -L$BREW_PREFIX/opt/libomp/lib -lomp \
      -Wl,-rpath,$BREW_PREFIX/opt/libomp/lib -o $tmpdir/conftest" >&2 || true
  rm -rf "$tmpdir"
  return 1
}

check_openmp_linux () {
  tmpdir="$(make_conftest_dir)" || return 1
  cat > "$tmpdir/conftest.cpp" <<'EOF'
#include <omp.h>
int main() {
  omp_set_num_threads(1);
  return omp_get_max_threads() >= 1 ? 0 : 1;
}
EOF
  if eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD -fopenmp \
      $tmpdir/conftest.cpp -o $tmpdir/conftest" \
      >/dev/null 2>&1; then
    rm -rf "$tmpdir"
    return 0
  fi
  echo "configure: Linux OpenMP link probe failed:" >&2
  eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD -fopenmp \
      $tmpdir/conftest.cpp -o $tmpdir/conftest" >&2 || true
  rm -rf "$tmpdir"
  return 1
}

detect_cppad_include_dir () {
  echo "configure: checking CppAD header locations..." >&2
  
  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists cppad 2>/dev/null; then
    result=$(pkg-config --cflags-only-I cppad 2>/dev/null \
      | sed -n 's/.*-I\([^[:space:]]*\).*/\1/p' \
      | head -n 1)
    echo "configure: pkg-config CppAD include path: ${result:-none}" >&2
    if [ -n "$result" ] && [ -f "$result/cppad/cppad.hpp" ]; then
      echo "configure: found CppAD headers via pkg-config at $result" >&2
      printf '%s' "$result"
      return 0
    fi
  fi

  for base in \
    "$BREW_PREFIX/opt/cppad/include" \
    /opt/homebrew/include \
    /usr/include \
    /usr/local/include \
    /opt/include
  do
    echo "configure: checking $base/cppad/cppad.hpp" >&2
    if [ -f "$base/cppad/cppad.hpp" ]; then
      echo "configure: found CppAD headers at $base" >&2
      printf '%s' "$base"
      return 0
    fi
  done

  echo "configure: ERROR: CppAD headers NOT found in any checked location" >&2
  return 1
}

detect_cppad_lib () {
  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists cppad 2>/dev/null; then
    pkg-config --libs cppad 2>/dev/null
    return 0
  fi

  for lib in \
    "$BREW_PREFIX/opt/cppad/lib/libcppad_lib.dylib" \
    /opt/homebrew/lib/libcppad_lib.dylib \
    /usr/local/lib/libcppad_lib.dylib \
    /opt/lib/libcppad_lib.so \
    /opt/homebrew/lib/libcppad_lib.so \
    /usr/local/lib/libcppad_lib.so
  do
    if [ -f "$lib" ]; then
      libdir="$(dirname "$lib")"
      if [ "$UNAME_S" = "Darwin" ]; then
        printf '%s' "$lib -Wl,-rpath,$libdir"
      else
        printf '%s' "-L$libdir -lcppad_lib"
      fi
      return 0
    fi
  done

  return 1
}

# ---- Brew discovery (macOS only) ----
if [ "$UNAME_S" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
  BREW_PREFIX="$(brew --prefix 2>/dev/null || true)"
  if [ -n "$BREW_PREFIX" ] && [ -d "$BREW_PREFIX/include" ]; then
    BREW_CPPFLAGS="-I$BREW_PREFIX/include"
  fi
fi

# Patch selected CppAD headers so the .so does not reference std::cout
# (R CMD check "compiled code" NOTE). Shadows live under src/include/cppad/
# and are removed by ./cleanup. See src/include/cppad/README.md.
generate_cppad_cout_shadows () {
  cppad_inc="$1"
  for rel in \
    local/var_op/atomic_op.hpp \
    core/ad_fun.hpp \
    core/fun_construct.hpp
  do
    src="$cppad_inc/cppad/$rel"
    dst="src/include/cppad/$rel"
    if [ ! -f "$src" ]; then
      echo "configure: ERROR: missing CppAD header needed for cout shadow: $src" >&2
      return 1
    fi
    mkdir -p "$(dirname "$dst")"
    {
      printf '%s\n' \
        "/* GENERATED by configure -- do not edit." \
        "   adlaplace R-safe shadow: std::cout -> ::adlaplace::cppad_trace_stream()" \
        "   Source: $src" \
        "   Revert: delete this file (also removed by ./cleanup). */" \
        "#include \"adlaplace/cppad_trace_stream.hpp\""
      sed -e 's/std::cout/::adlaplace::cppad_trace_stream()/g' "$src"
    } > "$dst"
  done
  echo "configure: generated CppAD std::cout shadow headers under src/include/cppad/" >&2
  return 0
}

CPPAD_INCLUDE_DIR="$(detect_cppad_include_dir || true)"
if [ -n "$CPPAD_INCLUDE_DIR" ]; then
  CPPAD_INCLUDE_CPPFLAGS="-I$CPPAD_INCLUDE_DIR"
  msg "configure: using CppAD headers from $CPPAD_INCLUDE_DIR"
  generate_cppad_cout_shadows "$CPPAD_INCLUDE_DIR" || exit 1
else
  echo "configure: ERROR: could not locate CppAD headers (cppad/cppad.hpp)" >&2
  exit 1
fi

CPPAD_LIBS="$(detect_cppad_lib || true)"
if [ -n "$CPPAD_LIBS" ]; then
  msg "configure: using CppAD library flags: $CPPAD_LIBS"
else
  msg "configure: no standalone CppAD library detected; relying on header-only CppAD"
fi

# ---- OpenMP (required; dynamic libomp on macOS; -fopenmp on Linux) ----
# Always report OpenMP status on stderr (visible in R CMD check install logs).
if [ "$UNAME_S" = "Darwin" ]; then
  if [ -n "$BREW_PREFIX" ] \
     && [ -f "$BREW_PREFIX/opt/libomp/include/omp.h" ] \
     && [ -f "$BREW_PREFIX/opt/libomp/lib/libomp.dylib" ]; then
    OPENMP_CPPFLAGS="-I$BREW_PREFIX/opt/libomp/include"
    OPENMP_CXXFLAGS="-Xpreprocessor -fopenmp"
    OPENMP_LIBS="-L$BREW_PREFIX/opt/libomp/lib -lomp -Wl,-rpath,$BREW_PREFIX/opt/libomp/lib"
    if check_openmp_darwin; then
      echo "configure: OpenMP enabled on macOS using dynamic Homebrew libomp" >&2
    else
      # Headers/dylib present: still enable flags. A failed probe can be
      # transient (TMPDIR/noexec, PATH) but omitting -lomp breaks dyn.load.
      echo "configure: WARNING: OpenMP link probe failed; enabling Homebrew libomp flags anyway" >&2
    fi
  else
    echo "configure: ERROR: OpenMP required but Homebrew libomp not found" >&2
    echo "configure: install with: brew install libomp" >&2
    exit 1
  fi
else
  if check_cxxflag "-fopenmp" && check_openmp_linux; then
    OPENMP_CXXFLAGS="-fopenmp"
    OPENMP_LIBS="-fopenmp"
    echo "configure: OpenMP enabled using -fopenmp" >&2
  else
    echo "configure: ERROR: OpenMP is required (SystemRequirements) but -fopenmp is not usable" >&2
    echo "configure: install a GCC/Clang toolchain with OpenMP (e.g. libgomp) and retry" >&2
    exit 1
  fi
fi

# Always true because we force it
msg "configure: CppAD ColPack: DISABLED (forced via CPPAD_HAS_COLPACK=0)"
