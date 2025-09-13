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

echo "=== Step 1: Generating k-core unweighted experiments ==="
./generate_experiments.sh
echo "Experiments generated successfully."
echo

echo "=== Step 2: Running k-core unweighted experiments ==="
./run_experiments.sh -n 1
echo "Experiments finished running."
echo

echo "=== Step 3: Extracting k-core unweighted results ==="
./extract_results.sh
echo "Results extracted successfully."
echo

echo "=== k-core unweighted pipeline completed successfully ==="
