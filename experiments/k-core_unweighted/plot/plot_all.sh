#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# ---------- CONFIG ----------
R_MAIN_SCRIPT="plot_unweighted.R"     # main R script to generate TikZ plots
TIKZ_MASTER="tikzplot.tex"      # wrapper LaTeX file
OUTPUT_NAME="kcore_unw"         # final PDF name
# -----------------------------

echo "[1/4] Converting CSVs via Python…"
python3 convert_all.py

echo "[2/4] Checking R deps (will auto-install if needed)…"
Rscript -e 'source("common.R"); cat("[deps] R package check complete.\n")'

echo "[3/4] Generating TikZ plots via R…"
Rscript "$R_MAIN_SCRIPT"

echo "[4/4] Compiling LaTeX wrapper to PDF…"
if ! command -v pdflatex >/dev/null 2>&1; then
  echo "[ERROR] pdflatex not found. Install TeX Live (or similar) to compile TikZ."
  exit 1
fi

# Run pdflatex in the same directory and output directly to OUTPUT_NAME
pdflatex -halt-on-error -interaction=batchmode -jobname="${OUTPUT_NAME%.pdf}" "$TIKZ_MASTER" >/dev/null

echo "✅ Done: $(realpath "$OUTPUT_NAME")"

rm kcore_unw.aux kcore_unw.log Rplots.pdf .cache.db