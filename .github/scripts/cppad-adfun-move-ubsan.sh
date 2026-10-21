#!/usr/bin/env bash
# UndefinedBehaviorSanitizer check for the vendored ADFun move constructor.
# Fetches unmodified CppAD 20260000.0, overlays RCppAD's fun_construct.hpp,
# and runs reproducers/cppad-adfun-move. A missing initializer aborts.
set -euo pipefail

root="${GITHUB_WORKSPACE:-}"
if [ -z "$root" ]; then
  root=$(cd "$(dirname "$0")/../.." && pwd)
fi
cxx="${CXX:-g++}"
bash "$root/reproducers/cppad-adfun-move/run.sh" --vendored --cxx "$cxx"
