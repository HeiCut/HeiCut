#!/bin/bash
###############################################################################
# perform_all_experiments.sh
#
# Runs the perform_<dataset>_experiments.sh script inside each experiment folder.
# The script cd's into each directory, executes the script, and returns back.
###############################################################################

# Base experiments folder (assuming this script is located in the same directory)
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# List of experiment subdirectories
SUBDIRS=(
  "medium_weighted"
  "medium_unweighted"
  "large_weighted"
  "large_unweighted"
  "k-core_weighted"
  "k-core_unweighted"
)

# Loop through each subdirectory and run the corresponding script
for dir in "${SUBDIRS[@]}"; do
  echo "======================================="
  echo "Entering $dir and running experiment..."
  echo "======================================="

  cd "$BASE_DIR/$dir" || { echo "Failed to enter $dir"; exit 1; }

  SCRIPT="perform_${dir}_experiments.sh"
  if [[ -x "$SCRIPT" ]]; then
    bash "$SCRIPT"
  else
    echo "Script $SCRIPT not found or not executable in $dir"
  fi

  # Return to base experiments folder
  cd "$BASE_DIR" || exit 1
done

echo "======================================="
echo "All experiments completed."
echo "======================================="
