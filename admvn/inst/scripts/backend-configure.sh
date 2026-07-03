# CppAD detection for admvn (header-only; optional libcppad_lib)

CPPAD_CPPFLAGS="-DCPPAD_HAS_COLPACK=0"
CPPAD_INCLUDE_CPPFLAGS=""
CPPAD_LIBS=""
BREW_CPPFLAGS=""
BREW_PREFIX=""

UNAME_S="$(uname -s 2>/dev/null || echo unknown)"

if [ "$UNAME_S" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
  BREW_PREFIX="$(brew --prefix 2>/dev/null || true)"
  if [ -n "$BREW_PREFIX" ] && [ -d "$BREW_PREFIX/include" ]; then
    BREW_CPPFLAGS="-I$BREW_PREFIX/include"
  fi
fi

detect_cppad_include_dir () {
  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists cppad 2>/dev/null; then
    result=$(pkg-config --cflags-only-I cppad 2>/dev/null \
      | sed -n 's/.*-I\([^[:space:]]*\).*/\1/p' \
      | head -n 1)
    if [ -n "$result" ] && [ -f "$result/cppad/cppad.hpp" ]; then
      printf '%s' "$result"
      return 0
    fi
  fi

  for base in \
    "$BREW_PREFIX/opt/cppad/include" \
    /opt/homebrew/include \
    /usr/local/include \
    /usr/include \
    /opt/include
  do
    if [ -f "$base/cppad/cppad.hpp" ]; then
      printf '%s' "$base"
      return 0
    fi
  done

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

CPPAD_INCLUDE_DIR="$(detect_cppad_include_dir || true)"
if [ -n "$CPPAD_INCLUDE_DIR" ]; then
  CPPAD_INCLUDE_CPPFLAGS="-I$CPPAD_INCLUDE_DIR"
else
  echo "configure: ERROR: could not locate CppAD headers (cppad/cppad.hpp)" >&2
  exit 1
fi

CPPAD_LIBS="$(detect_cppad_lib || true)"
