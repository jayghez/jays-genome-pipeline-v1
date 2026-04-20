#!/usr/bin/env bash
set -euo pipefail

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools is missing."
  echo "Install with: brew install bcftools htslib"
  echo "Or with conda: conda install -c bioconda bcftools htslib"
  exit 1
fi

echo "bcftools found: $(bcftools --version | head -n 1)"
