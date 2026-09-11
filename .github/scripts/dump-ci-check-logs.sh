#!/usr/bin/env bash
# Print R CMD check diagnostics from .ci-check/ (used when check-dir is
# ../.ci-check/... so r-lib/actions cannot upload or find "check/").
set -euo pipefail

root="${1:-.ci-check}"
if [ ! -d "$root" ]; then
  echo "No check directory: $root"
  exit 0
fi

echo "===== Listing $root ====="
find "$root" -maxdepth 4 \( -type d -o -name '00check.log' -o -name 'testthat.Rout*' \) -print | sort
echo

for rcheck in $(find "$root" -type d -name '*.Rcheck' | sort); do
  echo "===== $rcheck ====="
  for f in 00check.log 00install.out 00check-diff.log; do
    path="$rcheck/$f"
    if [ -f "$path" ]; then
      echo "----- $f -----"
      cat "$path"
      echo
    fi
  done
  if [ -d "$rcheck/tests" ]; then
    echo "----- tests/ listing -----"
    ls -la "$rcheck/tests" || true
    for rout in $(find "$rcheck/tests" -name 'testthat.Rout*' -type f | sort); do
      echo "----- $rout -----"
      cat "$rout"
      echo
    done
  fi
  echo
done
