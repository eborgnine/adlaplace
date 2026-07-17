# OpenMP detection for admvn (local copy of adlaplace/inst/scripts/backend-configure.sh).
# CppAD comes from LinkingTo: RCppAD (no system CppAD / -lcppad_lib).

OPENMP_CPPFLAGS=""
OPENMP_CXXFLAGS=""
OPENMP_LIBS=""
VISIBILITY_CXXFLAGS=""
BREW_PREFIX=""

UNAME_S="$(uname -s 2>/dev/null || echo unknown)"

# R CMD INSTALL/check set R_HOME; bare `R` is rejected there (R Installation §1.6).
# Use ${R_HOME:-} so `set -u` configure wrappers do not fail when unset.
if [ -z "${R_HOME:-}" ]; then
  R_HOME=`R RHOME`
fi
R_EXE="${R_HOME}/bin/R"
CXX_CMD="$("${R_EXE}" CMD config CXX 2>/dev/null || echo c++)"
CXXFLAGS_CMD="$("${R_EXE}" CMD config CXXFLAGS 2>/dev/null || true)"
CPPFLAGS_CMD="$("${R_EXE}" CMD config CPPFLAGS 2>/dev/null || true)"

make_conftest_dir () {
  mktemp -d "${TMPDIR:-/tmp}/admvn-conftest.XXXXXX"
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

# Do not add -I$BREW_PREFIX/include (would risk system CppAD over RCppAD).
if [ "$UNAME_S" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
  BREW_PREFIX="$(brew --prefix 2>/dev/null || true)"
fi

# ---- OpenMP (optional; SystemRequirements lists OpenMP as recommended) ----
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
      echo "configure: WARNING: OpenMP link probe failed; enabling Homebrew libomp flags anyway" >&2
    fi
  else
    echo "configure: OpenMP not enabled on macOS (Homebrew libomp not found)" >&2
  fi
else
  if check_cxxflag "-fopenmp" && check_openmp_linux; then
    OPENMP_CXXFLAGS="-fopenmp"
    OPENMP_LIBS="-fopenmp"
    echo "configure: OpenMP enabled using -fopenmp" >&2
  else
    echo "configure: OpenMP not enabled (compiler does not accept -fopenmp)" >&2
  fi
fi

# ---- Symbol visibility (keep CppAD statics local to this DSO) ----
# Only emit -fvisibility=hidden after a compiler probe (CRAN non-portable flags).
if check_cxxflag "-fvisibility=hidden"; then
  VISIBILITY_CXXFLAGS="-fvisibility=hidden"
  echo "configure: using -fvisibility=hidden" >&2
else
  echo "configure: -fvisibility=hidden not supported; leaving default visibility" >&2
fi
