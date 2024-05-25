#!/bin/bash

# Default R version
R_VERSION="4.1.2"

# Check if an R version is provided as an argument
if [ ! -z "$1" ]; then
  R_VERSION="$1"
fi

# Path to the Rscript binary
RSCRIPT_PATH="/opt/R/${R_VERSION}/bin/Rscript"

# Check if the Rscript binary exists
if [ ! -f "$RSCRIPT_PATH" ]; then
  echo "Rscript binary for R version ${R_VERSION} not found at ${RSCRIPT_PATH}"
  exit 1
fi

# Path to your R script
RSCRIPT_FILE="main.R"
shift;  # will remove first (which is the version) arg from the "$@"
# Execute the R script with the specified R version and forward any additional arguments
"$RSCRIPT_PATH" "$RSCRIPT_FILE" "$@"

