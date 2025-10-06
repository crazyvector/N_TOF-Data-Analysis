#!/bin/bash

# ====================================================================
# 2D Plot SCRIPT
# ====================================================================

# Check if 'jq' is installed and if a configuration file was provided
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it to run this script."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_2DPlot.sh config_2D.json"
    exit 1
fi

CONFIG_FILE="$1"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi

echo "=========================================================="
echo "Reading configuration from $CONFIG_FILE"
echo "=========================================================="

# Citim valorile din JSON folosind jq
ROOT_FOLDER=$(jq -r '.rootFile' "$CONFIG_FILE")
CALIB_FILE=$(jq -r '.calibrationFile' "$CONFIG_FILE")
DETECTOR=$(jq -r '.detectorName' "$CONFIG_FILE")
DETECTOR_COUNT=$(jq -r '.detectorCount' "$CONFIG_FILE")
SOURCE=$(jq -r '.sourceName' "$CONFIG_FILE")
FIT_TYPE=$(jq -r '.fitType' "$CONFIG_FILE")
PS_INT_THRESHOLD=$(jq -r '.ps_int_threshold' "$CONFIG_FILE")
DISTANCE=$(jq -r '.distance' "$CONFIG_FILE")

# X Axis
TMIN=$(jq -r '.xAxis.tMin' "$CONFIG_FILE")
TMAX=$(jq -r '.xAxis.tMax' "$CONFIG_FILE")
NDEC=$(jq -r '.xAxis.ndec' "$CONFIG_FILE")
NBPDEC=$(jq -r '.xAxis.nbpdec' "$CONFIG_FILE")
NBINSX=$(jq -r '.xAxis.nbinsX' "$CONFIG_FILE")
STEP=$(jq -r '.xAxis.step' "$CONFIG_FILE")

# Y Axis
YMIN=$(jq -r '.yAxis.yMin' "$CONFIG_FILE")
YMAX=$(jq -r '.yAxis.yMax' "$CONFIG_FILE")
NBINY=$(jq -r '.yAxis.n_bin_y' "$CONFIG_FILE")
BINWIDTH=$(jq -r '.yAxis.binWidth' "$CONFIG_FILE")

# Citim lista de runs
RUNS_ARRAY=($(jq -r '.runs[]' "$CONFIG_FILE"))

# Construim lista ROOT file paths
ROOT_FILES=()
for RUN in "${RUNS_ARRAY[@]}"; do
    ROOT_FILES+=("$ROOT_FOLDER/run$RUN.root")
done

# Afisam pentru verificare
echo "ROOT folder: $ROOT_FOLDER"
echo "Calibration file: $CALIB_FILE"
echo "Detector: $DETECTOR ($DETECTOR_COUNT)"
echo "Source: $SOURCE"
echo "Fit type: $FIT_TYPE"
echo "Runs: ${RUNS_ARRAY[*]}"
echo "X Axis: TMIN=$TMIN, TMAX=$TMAX, NDEC=$NDEC, NBPDEC=$NBPDEC, NBINSX=$NBINSX, STEP=$STEP"
echo "Y Axis: YMIN=$YMIN, YMAX=$YMAX, NBINY=$NBINY, BINWIDTH=$BINWIDTH"

# Construim comanda pentru ROOT
ROOT_CMD=".L Plot2D.C+
plot2D({$(printf '"%s",' "${RUNS_ARRAY[@]}" | sed 's/,$//')}, \"$DETECTOR\", $DETECTOR_COUNT, \"$SOURCE\", \"$FIT_TYPE\", $TMIN, $TMAX, $NDEC, $NBPDEC, $NBINSX, $STEP, $YMIN, $YMAX, $NBINY, $BINWIDTH, \"$CALIB_FILE\", \"$ROOT_FOLDER\", $PS_INT_THRESHOLD, $DISTANCE);"

echo "Executing ROOT command:"
echo "$ROOT_CMD"

# Apelam ROOT
echo "$ROOT_CMD" | root -l -b
