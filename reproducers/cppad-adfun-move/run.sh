#!/usr/bin/env bash
# Compile and run adfun_move.cpp under UndefinedBehaviorSanitizer.
#
# Against unmodified CppAD 20260000.0 (downloads the tag, then adds the
# generated configure.hpp that the tag archive does not ship):
#   ./run.sh --fetch 20260000.0 --cxx clang++
#
# Same tree with only the ADFun move constructor initialized before swap:
#   ./run.sh --fetch 20260000.0 --patch-move --cxx clang++
#
# Upstream tag plus the vendored, already patched fun_construct.hpp:
#   ./run.sh --vendored --cxx clang++
#
# CPPAD_INCLUDE and CXX are accepted in place of --include and --cxx.
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
include="${CPPAD_INCLUDE:-}"
cxx="${CXX:-c++}"
fetch=""
patch_move=0
vendored=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    --include)
      include="${2:-}"
      shift 2
      ;;
    --cxx)
      cxx="${2:-}"
      shift 2
      ;;
    --fetch)
      fetch="${2:-}"
      shift 2
      ;;
    --patch-move)
      patch_move=1
      shift
      ;;
    --vendored)
      vendored=1
      shift
      ;;
    *)
      echo "unknown argument: $1" >&2
      echo "usage: $0 (--include DIR | --fetch VERSION | --vendored) [--patch-move] [--cxx COMPILER]" >&2
      exit 2
      ;;
  esac
done

if [ "$vendored" -eq 1 ] && [ "$patch_move" -eq 1 ]; then
  echo "pass only one of --vendored and --patch-move" >&2
  exit 2
fi
if [ "$vendored" -eq 1 ]; then
  fetch="${fetch:-20260000.0}"
fi
if [ -n "$fetch" ] && [ -n "$include" ]; then
  echo "pass only one of --include and --fetch/--vendored" >&2
  exit 2
fi

work="$(mktemp -d "${TMPDIR:-/tmp}/adfun-move.XXXXXX")"
trap 'rm -rf "$work"' EXIT

if [ -n "$fetch" ]; then
  cache="${TMPDIR:-/tmp}/cppad-${fetch}"
  if [ ! -f "$cache/include/cppad/cppad.hpp" ]; then
    mkdir -p "$cache"
    url="https://github.com/coin-or/CppAD/archive/refs/tags/${fetch}.tar.gz"
    echo "fetch: $url"
    curl -fsSL "$url" | tar -xz -C "$cache" --strip-components=1
  fi
  if [ ! -f "$cache/include/cppad/core/fun_construct.hpp" ]; then
    echo "downloaded CppAD ${fetch} has no include/cppad/core/fun_construct.hpp" >&2
    exit 2
  fi
  # The release tag ships configure.hpp.in. Overlay a generated configure.hpp
  # without modifying fun_construct.hpp.
  stage="$work/upstream-include"
  mkdir -p "$stage"
  cp -a "$cache/include/." "$stage/"
  cp "$root/configure.hpp" "$stage/cppad/configure.hpp"
  include="$stage"
fi

if [ "$vendored" -eq 1 ]; then
  repo="$(cd "$root/../.." && pwd)"
  vendored_ctor="$repo/RCppAD/inst/include/cppad/core/fun_construct.hpp"
  trace_header="$repo/RCppAD/inst/include/RCppAD/cppad_trace_stream.hpp"
  if [ ! -f "$vendored_ctor" ] || [ ! -f "$trace_header" ]; then
    echo "vendored CppAD headers not found under $repo/RCppAD" >&2
    exit 2
  fi
  cp "$vendored_ctor" "$include/cppad/core/fun_construct.hpp"
  mkdir -p "$include/RCppAD"
  cp "$trace_header" "$include/RCppAD/cppad_trace_stream.hpp"
fi

if [ "$patch_move" -eq 1 ]; then
  case "$include" in
    "$work"/*) ;;
    *)
      patched="$work/patched-include"
      mkdir -p "$patched"
      cp -a "$include/." "$patched/"
      include="$patched"
      ;;
  esac
  ctor="$include/cppad/core/fun_construct.hpp"
  tmp="$ctor.patched"
  awk '
    /none of the default constructor values matter to the destructor/ {
      print "// Initialize before swap. std::swap loads the destination bools."
      next
    }
    /::ADFun\(ADFun&& f\)/ {
      print
      print ": ADFun<Base, RecBase>()"
      next
    }
    { print }
  ' "$ctor" > "$tmp"
  mv "$tmp" "$ctor"
  if ! grep -A1 '::ADFun(ADFun&& f)' "$ctor" | grep -q ': ADFun<Base, RecBase>()'; then
    echo "failed to patch $ctor" >&2
    exit 2
  fi
fi

if [ -z "$include" ]; then
  echo "pass --include, --fetch VERSION, or set CPPAD_INCLUDE" >&2
  exit 2
fi
if [ ! -f "$include/cppad/cppad.hpp" ]; then
  echo "missing $include/cppad/cppad.hpp" >&2
  exit 2
fi
if [ ! -f "$include/cppad/configure.hpp" ]; then
  echo "missing $include/cppad/configure.hpp (use --fetch for a tag archive)" >&2
  exit 2
fi

cmd=(
  "$cxx" -std=c++17 -O1 -g
  -fsanitize=undefined -fno-sanitize-recover=undefined
  -I"$include"
  "$root/adfun_move.cpp" -o "$work/adfun_move"
)
printf 'compile:'
printf ' %q' "${cmd[@]}"
printf '\n'
"${cmd[@]}"
printf 'run: %q\n' "$work/adfun_move"
"$work/adfun_move"
