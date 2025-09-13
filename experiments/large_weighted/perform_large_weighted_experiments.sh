#!/bin/bash
###############################################################################
# master_pipeline.sh
# #
# Runs the full experiment pipeline:
#   1) generate_experiments.sh
#   2) run_experiments.sh
#   3) extract_results.sh
# #
###############################################################################

set -euo pipefail

echo "=== Step 1: Generating large weighted experiments ==="
./generate_experiments.sh
echo "Experiments generated successfully."
echo

echo "=== Step 2: Running large weighted experiments ==="
./run_experiments.sh -n 1
echo "Experiments finished running."
echo

echo "=== Step 3: Extracting large weighted results ==="
./extract_results.sh
echo "Results extracted successfully."
echo

echo "=== large weighted pipeline completed successfully ==="
