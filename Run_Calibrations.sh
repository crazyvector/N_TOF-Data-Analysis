#!/bin/bash

# ====================================================================
# ENERGY CALIBRATION SCRIPT
# This script automates the calibration process for multiple
# runs and detectors, handling different sources and fit types.
# ====================================================================

# --- Step 1: Check dependencies ---
# Verify if 'jq' is installed, which is required to read JSON files.
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it to run this script."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

# --- Step 2: Check for JSON config argument ---
# The script expects a JSON configuration file as its first argument.
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_Calibration.sh config.json"
    exit 1
fi

CONFIG_FILE="$1"

# --- Step 3: Validate JSON file existence ---
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi

echo "=========================================================="
echo "Reading configuration from $CONFIG_FILE"
echo "=========================================================="

# --- Step 4: Clear previous parameters ---
# Ensure calibration_parameters.txt is empty before writing new values.
> calibration_parameters.txt

# ===============================================
# 1. Read global and advanced parameters from JSON
# ===============================================

# Input and output folders
INPUT_FOLDER=$(jq -r '.inputFolder' "$CONFIG_FILE")
OUTPUT_FOLDER=$(jq -r '.outputFolder' "$CONFIG_FILE")
FIT_TYPE=$(jq -r '.fitType' "$CONFIG_FILE")

# Histogram parameters
REBIN_FACTOR=$(jq -r '.rebinFactor' "$CONFIG_FILE")
XMIN_HIST=$(jq -r '.xminHist' "$CONFIG_FILE")
XMAX_HIST=$(jq -r '.xmaxHist' "$CONFIG_FILE")
BIN_WIDTH=$(jq -r '.binWidth' "$CONFIG_FILE")

# Ensure paths end with a slash for consistency
if [[ "$INPUT_FOLDER" != */ ]]; then
    INPUT_FOLDER="$INPUT_FOLDER/"
fi
if [[ "$OUTPUT_FOLDER" != */ ]]; then
    OUTPUT_FOLDER="$OUTPUT_FOLDER/"
fi

# ===============================================
# 2. Iterate through sources in JSON
# ===============================================

echo "=========================================================="
echo "Processing sources..."

# Loop over each source defined in the JSON
for i in $(jq -r '.sources | keys[]' "$CONFIG_FILE"); do
    SOURCE_DATA=$(jq -r ".sources[$i]" "$CONFIG_FILE")

    # --- Extract source-specific information ---
    n_runs=$(echo "$SOURCE_DATA" | jq -r '.n_runs')                          # Number of runs
    run_names_str=$(echo "$SOURCE_DATA" | jq -r '.run_names | join(" ")')     # Run names as string
    detector=$(echo "$SOURCE_DATA" | jq -r '.detector')                       # Detector tree name
    ndet=$(echo "$SOURCE_DATA" | jq -r '.ndet')                               # Number of detectors
    source=$(echo "$SOURCE_DATA" | jq -r '.source')                           # Source name
    ids_str=$(echo "$SOURCE_DATA" | jq -r '.ids | join(" ")')                 # Detector IDs for fitting

    # --- Extract fit parameters specific to each source ---
    energies_str=$(echo "$SOURCE_DATA" | jq -r '.energies | join(" ")')       # Known gamma energies
    FIT_RATIO=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.fitRatio')
    FIT_SIGMA=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.fitSigma')
    PLOT_MIN_RANGE=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.plotMinRange')
    PLOT_MAX_RANGE=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.plotMaxRange')
    C1=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.c1')
    C2=$(echo "$SOURCE_DATA" | jq -r '.fitParameters.c2')

    # --- Use global histogram parameters ---
    rebin_factor="$REBIN_FACTOR"
    ratio="$FIT_RATIO"
    sigma="$FIT_SIGMA"
    range_xmin="$PLOT_MIN_RANGE"
    range_xmax="$PLOT_MAX_RANGE"
    c1="$C1"
    c2="$C2"
    xmin_hist="$XMIN_HIST"
    xmax_hist="$XMAX_HIST"
    bin_width="$BIN_WIDTH"

    # --- Convert space-separated strings to Bash arrays ---
    read -ra run_names <<< "$run_names_str"
    read -ra ids <<< "$ids_str"
    read -ra energies <<< "$energies_str"

    # --- Print source details for verification ---
    echo "=========================================================="
    echo "Source details:"
    echo "  Runs: ${run_names[*]}"
    echo "  Detector: $detector"
    echo "  Number of detectors: $ndet"
    echo "  Source name: $source"
    echo "  IDs for fitting: ${ids[*]}"
    echo "  Energies for fitting: ${energies[*]}"
    echo "  Rebin factor: $rebin_factor"
    echo "  Fit parameters: Ratio=$ratio, Sigma=$sigma, Range=[$range_xmin, $range_xmax], Constants=[$c1, $c2]"

    # ===============================================
    # Calibration 1: Sum or individual histograms
    # ===============================================
    echo "==> Running Calibration 1 for all selected runs"

    # Determine the output histogram file name
    histfile="${OUTPUT_FOLDER}Rebin_${source}_summed_histograms.root"
    if [[ "$n_runs" -eq 1 ]]; then
        histfile="${OUTPUT_FOLDER}histograma_${run_names[0]}_${detector}_${source}.root"
    fi

    # Convert Bash array to C++ vector string
    run_cpp=$(printf ', "%s"' "${run_names[@]}")
    run_cpp="{${run_cpp:2}}"

    if [[ -f "$histfile" ]]; then
        echo "⚠️  Histogram $histfile already exists. Skipping Calibration 1."
    else
        # Build ROOT command for Calibration 1
        ROOT_CMD='.L 1_CalibrationsD.C
        calibrations('"$run_cpp"', "'$detector'", '"$ndet"', "'$source'", "'$INPUT_FOLDER'", "'$OUTPUT_FOLDER'",'"$rebin_factor"', '"$xmin_hist"', '"$xmax_hist"', '"$bin_width"');
        '
        echo "Executing ROOT command for Calibration 1:"
        echo "$ROOT_CMD"
        echo "$ROOT_CMD" | root -l -b
    fi

    # ===============================================
    # Calibration 2: Fit Peaks for energy calibration
    # ===============================================
    # Build C++ vector string for detector IDs
    ids_cpp=$(printf ",%s" "${ids[@]}")
    ids_cpp="{${ids_cpp:1}}"

    # Build C++ vector string for energies
    energies_cpp=$(printf ",%s" "${energies[@]}")
    energies_cpp="{${energies_cpp:1}}"

    echo "=========================================================="
    echo "==> Running Calibration 2 for source $source with IDs: ${ids[*]}"

    ROOT_CMD='.L 2_CalibrationsD.C
    calibrations("'"$histfile"'", "'"$detector"'", '"$ndet"', "'"$source"'", '"$ids_cpp"', '"$energies_cpp"', '"$ratio"', '"$sigma"', '"$range_xmin"', '"$range_xmax"', '"$c1"', '"$c2"', "'"$OUTPUT_FOLDER"'");
    '
    echo "Executing ROOT command for Calibration 2:"
    echo "$ROOT_CMD"
    echo "$ROOT_CMD" | root -l -b

done

# ===============================================
# Calibration 3: Fit type selection for all detectors
# ===============================================
echo "=========================================================="
echo "==> Fit for detector(s) ${ids[*]} with fit type '$FIT_TYPE'"
filename="${OUTPUT_FOLDER}Edge_Compton_Detector="$detector".txt"

ROOT_CMD='.L 3_CalibrationsD.C
fit("'$filename'", "'$FIT_TYPE'", "'$OUTPUT_FOLDER'");
'
echo "Executing ROOT command for Calibration 3:"
echo "$ROOT_CMD"
echo "$ROOT_CMD" | root -l -b

# ===============================================
# Calibration 4: Convert Channel -> Energy for each source
# ===============================================
echo "=========================================================="
echo "==> Running Channel to Energy for all sources"
# Re-iterate through sources to handle Channel->Energy conversion
for i in $(jq -r '.sources | keys[]' "$CONFIG_FILE"); do
    SOURCE_DATA=$(jq -r ".sources[$i]" "$CONFIG_FILE")
    run_names_str=$(echo "$SOURCE_DATA" | jq -r '.run_names | join(" ")')
    read -ra run_names <<< "$run_names_str"
    source=$(echo "$SOURCE_DATA" | jq -r '.source')

    run_cpp=$(printf ', "%s"' "${run_names[@]}")
    run_cpp="{${run_cpp:2}}"

    echo "Source: $source"
    echo "Runs:   ${run_names[*]}"

    ROOT_CMD=".L 4_CalibrationsD.C
calibrations_energy($run_cpp, \"$detector\", $ndet, \"$source\", \"$FIT_TYPE\", \"$INPUT_FOLDER\", \"$OUTPUT_FOLDER\", $xmin_hist, $xmax_hist, $bin_width);"

    echo "Executing ROOT command for Calibration 4:"
    echo "$ROOT_CMD"
    echo "$ROOT_CMD" | root -l -b
done

echo "=========================================================="
echo "✅ All scripts executed."
echo "=========================================================="
