#!/usr/bin/env bash
# Compare vendored CppAD (RCppAD) to the latest coin-or/CppAD GitHub release.
# Soft notice only — never fails the job.
set -euo pipefail

ROOT="$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)"
CONF="$ROOT/RCppAD/inst/include/cppad/configure.hpp"

if [[ ! -f "$CONF" ]]; then
  echo "check-cppad-upstream: missing $CONF" >&2
  exit 0
fi

local_ver="$(sed -n 's/.*"cppad-\([^"]*\)".*/\1/p' "$CONF" | head -n 1)"
if [[ -z "$local_ver" ]]; then
  echo "check-cppad-upstream: could not parse CPPAD_PACKAGE_STRING" >&2
  exit 0
fi

api_url="https://api.github.com/repos/coin-or/CppAD/releases/latest"
curl_args=(
  -fsSL
  -H "Accept: application/vnd.github+json"
  -H "X-GitHub-Api-Version: 2022-11-28"
)
token="${GITHUB_PAT:-${GITHUB_TOKEN:-}}"
if [[ -n "$token" ]]; then
  curl_args+=(-H "Authorization: Bearer $token")
fi

json="$(curl "${curl_args[@]}" "$api_url" 2>/dev/null)" || {
  echo "check-cppad-upstream: could not fetch $api_url (network?); skipping" >&2
  exit 0
}

latest_ver="$(printf '%s\n' "$json" | sed -n 's/.*"tag_name"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' | head -n 1)"
if [[ -z "$latest_ver" ]]; then
  echo "check-cppad-upstream: no tag_name in releases/latest; skipping" >&2
  exit 0
fi

echo "Vendored CppAD: $local_ver"
echo "Latest coin-or/CppAD release: $latest_ver"

if [[ "$local_ver" == "$latest_ver" ]]; then
  echo "RCppAD CppAD headers are up to date with the latest GitHub release."
  exit 0
fi

# Numeric-ish compare for CppAD tags like 20260000.0
newer="$(printf '%s\n%s\n' "$local_ver" "$latest_ver" | sort -V | tail -n 1)"
if [[ "$newer" == "$latest_ver" && "$local_ver" != "$latest_ver" ]]; then
  msg="Upstream CppAD $latest_ver is newer than vendored $local_ver. Refresh with: cd RCppAD && ./tools/update-cppad.sh $latest_ver"
  echo "$msg"
  # GitHub Actions annotation (ignored outside GHA)
  echo "::notice title=CppAD upstream available::$msg"
else
  echo "Vendored CppAD $local_ver is ahead of (or incomparable to) release $latest_ver."
fi

exit 0
