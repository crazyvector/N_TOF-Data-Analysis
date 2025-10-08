#!/bin/bash
# ====================================================================
# 2D Plot SCRIPT
# This script reads a JSON configuration file and runs a ROOT macro
# to create 2D histograms (Neutron Energy vs Gamma Energy).
# ====================================================================

# --------------------------------------------------------------------
# Check if 'jq' is installed
# 'jq' is a command-line JSON processor used to parse the config file.
# --------------------------------------------------------------------
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it to run this script."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1   # exit with error code if jq is missing
fi

# --------------------------------------------------------------------
# Check if a configuration file argument is provided
# --------------------------------------------------------------------
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_2DPlot.sh config_2D.json"
    exit 1
fi

# Assign the first argument as the configuration file
CONFIG_FILE="$1"

# Check if the file actually exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi

# Print header to indicate reading config
echo "=========================================================="
echo "Reading configuration from $CONFIG_FILE"
echo "=========================================================="

# --------------------------------------------------------------------
# Read values from JSON configuration using 'jq'
# The '-r' flag ensures raw output (no quotes around strings)
# --------------------------------------------------------------------
ROOT_FOLDER=$(jq -r '.rootFile' "$CONFIG_FILE")        # folder containing ROOT files
CALIB_FILE=$(jq -r '.calibrationFile' "$CONFIG_FILE")  # calibration file path
DETECTOR=$(jq -r '.detectorName' "$CONFIG_FILE")       # detector tree name
DETECTOR_COUNT=$(jq -r '.detectorCount' "$CONFIG_FILE") # number of detectors
SOURCE=$(jq -r '.sourceName' "$CONFIG_FILE")           # source name
FIT_TYPE=$(jq -r '.fitType' "$CONFIG_FILE")            # type of calibration fit
PS_INT_THRESHOLD=$(jq -r '.ps_int_threshold' "$CONFIG_FILE") # pulse intensity threshold
CONDITION=$(jq -r '.condition' "$CONFIG_FILE")         # condition for pulse intensity filtering
DISTANCE=$(jq -r '.distance' "$CONFIG_FILE")           # distance for neutron energy calculation

# X Axis parameters
TMIN=$(jq -r '.xAxis.tMin' "$CONFIG_FILE")            # min log10 value
TMAX=$(jq -r '.xAxis.tMax' "$CONFIG_FILE")            # max log10 value
NDEC=$(jq -r '.xAxis.ndec' "$CONFIG_FILE")            # number of decades (for binning)
NBPDEC=$(jq -r '.xAxis.nbpdec' "$CONFIG_FILE")        # number of bins per decade
NBINSX=$(jq -r '.xAxis.nbinsX' "$CONFIG_FILE")        # total number of X bins
STEP=$(jq -r '.xAxis.step' "$CONFIG_FILE")            # logarithmic step size

# Y Axis parameters
YMIN=$(jq -r '.yAxis.yMin' "$CONFIG_FILE")            # minimum Y value
YMAX=$(jq -r '.yAxis.yMax' "$CONFIG_FILE")            # maximum Y value
NBINY=$(jq -r '.yAxis.n_bin_y' "$CONFIG_FILE")       # number of Y bins
BINWIDTH=$(jq -r '.yAxis.binWidth' "$CONFIG_FILE")    # Y bin width

# --------------------------------------------------------------------
# Read the list of runs from the JSON file
# The runs array is converted into a bash array
# --------------------------------------------------------------------
RUNS_ARRAY=($(jq -r '.runs[]' "$CONFIG_FILE"))

# --------------------------------------------------------------------
# Build full ROOT file paths using the ROOT_FOLDER and run numbers
# --------------------------------------------------------------------
ROOT_FILES=()
for RUN in "${RUNS_ARRAY[@]}"; do
    ROOT_FILES+=("$ROOT_FOLDER/run$RUN.root")
done

# --------------------------------------------------------------------
# Display configuration values for verification
# --------------------------------------------------------------------
echo "ROOT folder: $ROOT_FOLDER"
echo "Calibration file: $CALIB_FILE"
echo "Detector: $DETECTOR ($DETECTOR_COUNT)"
echo "Source: $SOURCE"
echo "Fit type: $FIT_TYPE"
echo "Runs: ${RUNS_ARRAY[*]}"
echo "X Axis: TMIN=$TMIN, TMAX=$TMAX, NDEC=$NDEC, NBPDEC=$NBPDEC, NBINSX=$NBINSX, STEP=$STEP"
echo "Y Axis: YMIN=$YMIN, YMAX=$YMAX, NBINY=$NBINY, BINWIDTH=$BINWIDTH"

# --------------------------------------------------------------------
# Construct the ROOT command to run the Plot2D.C macro
# 1. Load the macro with ".L Plot2D.C+"
#    the '+' compiles the macro for faster execution.
# 2. Call the plot2D() function with all parameters extracted from JSON.
# 3. The runs array is formatted into a ROOT-style { "run1", "run2", ... }.
# --------------------------------------------------------------------
ROOT_CMD=".L Plot2D.C
plot2D({$(printf '"%s",' "${RUNS_ARRAY[@]}" | sed 's/,$//')}, \"$DETECTOR\", $DETECTOR_COUNT, \"$SOURCE\", \"$FIT_TYPE\", $TMIN, $TMAX, $NDEC, $NBPDEC, $NBINSX, $STEP, $YMIN, $YMAX, $NBINY, $BINWIDTH, \"$CALIB_FILE\", \"$ROOT_FOLDER\", $PS_INT_THRESHOLD, $DISTANCE, \"$CONDITION\");"

# Print the ROOT command before execution
echo "Executing ROOT command:"
echo "$ROOT_CMD"

# --------------------------------------------------------------------
# Execute the ROOT command in batch mode (-b) and without splash screen (-l)
# --------------------------------------------------------------------
echo "$ROOT_CMD" | root -l -b
