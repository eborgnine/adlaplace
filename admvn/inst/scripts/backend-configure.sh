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

check_cxxflag () {
  flag="$1"
  tmpdir="config.$$"
  mkdir "$tmpdir"
  cat > "$tmpdir/conftest.cpp" <<'EOF'
int main() { return 0; }
EOF
  if eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD $flag -c $tmpdir/conftest.cpp -o $tmpdir/conftest.o" \
       >/dev/null 2>&1; then
    rm -rf "$tmpdir"
    return 0
  fi
  rm -rf "$tmpdir"
  return 1
}

check_openmp_darwin () {
  tmpdir="config.$$"
  mkdir "$tmpdir"
  cat > "$tmpdir/conftest.cpp" <<'EOF'
#include <omp.h>
int main() {
  omp_set_num_threads(1);
  return omp_get_max_threads() >= 1 ? 0 : 1;
}
EOF
  if eval "$CXX_CMD $CPPFLAGS_CMD $CXXFLAGS_CMD -I$BREW_PREFIX/opt/libomp/include -Xpreprocessor -fopenmp \
      $tmpdir/conftest.cpp $BREW_PREFIX/opt/libomp/lib/libomp.a -o $tmpdir/conftest" \
      >/dev/null 2>&1; then
    rm -rf "$tmpdir"
    return 0
  fi
  rm -rf "$tmpdir"
  return 1
}

check_openmp_linux () {
  tmpdir="config.$$"
  mkdir "$tmpdir"
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

CPPAD_INCLUDE_DIR="$(detect_cppad_include_dir || true)"
if [ -n "$CPPAD_INCLUDE_DIR" ]; then
  CPPAD_INCLUDE_CPPFLAGS="-I$CPPAD_INCLUDE_DIR"
  msg "configure: using CppAD headers from $CPPAD_INCLUDE_DIR"
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

# ---- Optional OpenMP (static libomp.a on macOS; dynamic -fopenmp on Linux) ----
if [ "$UNAME_S" = "Darwin" ]; then
  if [ -n "$BREW_PREFIX" ] \
     && [ -f "$BREW_PREFIX/opt/libomp/include/omp.h" ] \
     && [ -f "$BREW_PREFIX/opt/libomp/lib/libomp.a" ]; then
    if check_openmp_darwin; then
      OPENMP_CPPFLAGS="-I$BREW_PREFIX/opt/libomp/include"
      OPENMP_CXXFLAGS="-Xpreprocessor -fopenmp"
      OPENMP_LIBS="$BREW_PREFIX/opt/libomp/lib/libomp.a"
      msg "configure: OpenMP enabled on macOS using static Homebrew libomp"
    else
      msg "configure: OpenMP not enabled (compiler does not accept -Xpreprocessor -fopenmp)"
    fi
  else
    msg "configure: OpenMP not enabled on macOS (Homebrew libomp not found)"
  fi
else
  if check_cxxflag "-fopenmp" && check_openmp_linux; then
    OPENMP_CXXFLAGS="-fopenmp"
    OPENMP_LIBS="-fopenmp"
    msg "configure: OpenMP enabled using -fopenmp"
  else
    msg "configure: OpenMP not enabled (compiler does not accept -fopenmp)"
  fi
fi

# Always true because we force it
msg "configure: CppAD ColPack: DISABLED (forced via CPPAD_HAS_COLPACK=0)"
