#!/usr/bin/env bash
# Stage built vignette HTML into _site/<pkg>/ for GitHub Pages.
# Prefers vignettes from R CMD check trees (install uses --no-build-vignettes).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SITE="${ROOT}/_site"
rm -rf "$SITE"
mkdir -p "$SITE"

# Packages with vignettes that are checked in CI (plain check-dir).
CHECKED_PKGS=(adlaplace adlaplaceExample adlaplaceFem)

copy_doc_dir() {
  local pkg=$1
  local doc=$2
  local dest="${SITE}/${pkg}"
  if [ ! -d "$doc" ]; then
    return 1
  fi
  # Need at least one HTML vignette
  if ! compgen -G "${doc}/*.html" > /dev/null; then
    return 1
  fi
  mkdir -p "$dest"
  # HTML plus knitr companions (css, figures, etc.); skip R sources / rebuild inputs
  find "$doc" -maxdepth 1 -type f ! -name '*.R' ! -name '*.Rnw' ! -name '*.Rmd' \
    -exec cp {} "$dest/" \;
  # Figure / asset subdirs if present
  find "$doc" -mindepth 1 -maxdepth 1 -type d -exec cp -R {} "$dest/" \;
  echo "staged ${pkg} from ${doc}"
  return 0
}

stage_from_check() {
  local pkg=$1
  # CI writes check trees under repo-root .ci-check/ (never inside package dirs).
  local candidates=(
    "${ROOT}/.ci-check/${pkg}-plain/${pkg}.Rcheck/${pkg}/doc"
    "${ROOT}/.ci-check/${pkg}-as-cran/${pkg}.Rcheck/${pkg}/doc"
  )
  local doc
  for doc in "${candidates[@]}"; do
    if copy_doc_dir "$pkg" "$doc"; then
      return 0
    fi
  done
  echo "warning: no check vignettes found for ${pkg}" >&2
  return 1
}

staged=()
for pkg in "${CHECKED_PKGS[@]}"; do
  if stage_from_check "$pkg"; then
    staged+=("$pkg")
  fi
done

if [ "${#staged[@]}" -eq 0 ]; then
  echo "error: no vignette HTML staged" >&2
  exit 1
fi

# Root index listing packages and HTML vignettes
{
  cat <<'EOF'
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>adlaplace vignettes</title>
  <style>
    body { font-family: system-ui, sans-serif; max-width: 48rem; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; }
    h1 { font-size: 1.5rem; }
    h2 { font-size: 1.15rem; margin-top: 1.75rem; }
    ul { padding-left: 1.2rem; }
    a { color: #0b57d0; }
  </style>
</head>
<body>
  <h1>adlaplace vignettes</h1>
  <p>HTML vignettes built by CI from this monorepo.</p>
EOF

  for pkg in "${staged[@]}"; do
    echo "  <h2>${pkg}</h2>"
    echo "  <ul>"
    # Prefer listing *.html that look like vignettes (skip index.html if any)
    while IFS= read -r html; do
      base=$(basename "$html")
      [ "$base" = "index.html" ] && continue
      title=$base
      echo "    <li><a href=\"${pkg}/${base}\">${title}</a></li>"
    done < <(find "${SITE}/${pkg}" -maxdepth 1 -type f -name '*.html' | sort)
    echo "  </ul>"
  done

  cat <<'EOF'
</body>
</html>
EOF
} > "${SITE}/index.html"

echo "staged site at ${SITE} (packages: ${staged[*]})"
