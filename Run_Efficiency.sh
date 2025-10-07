#!/bin/bash

# ====================================================================
# EFFICIENCY CALCULATION SCRIPT
# This script automates the calculation of detector efficiency for
# multiple sources, reading parameters from a JSON configuration file
# and executing ROOT macros (PeakIntegrator.C and FitEfficiency.C).
# ====================================================================

# ----------------------------------------------------
# A. Initial checks
# ----------------------------------------------------
# Check if 'jq' is installed
# 'jq' is used to parse JSON configuration files
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

# Check if a JSON configuration file was provided as the first argument
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_Efficiency.sh config_eff.json"
    exit 1
fi

CONFIG_FILE="$1"

# Verify that the JSON configuration file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi

echo "=========================================================="
echo "Reading configuration from $CONFIG_FILE"
echo "=========================================================="

# ----------------------------------------------------
# B. Global parameters
# ----------------------------------------------------
# Read the fit model choice (usually 4 or 6) from JSON
FIT_MODEL=$(jq -r '.modelChoice' "$CONFIG_FILE")

# Read the energy interpolation type
ENERGY_INTERPOLATION=$(jq -r '.energyInterpolation' "$CONFIG_FILE")

# Read the detector IDs as a comma-separated string
IDS=$(jq -r '.ids | join(",")' "$CONFIG_FILE")

# Check if the fit model is valid; default to 4 if not
if [ "$FIT_MODEL" != "4" ] && [ "$FIT_MODEL" != "6" ]; then
    echo "⚠️ Fit model ('modelChoice') is not 4 or 6. Using 4 instead."
    FIT_MODEL=4
fi

# ----------------------------------------------------
# C. Iterate through efficiency sources
# ----------------------------------------------------
# Count the number of sources defined in JSON
NUM_SOURCES=$(jq -r '.sources | length' "$CONFIG_FILE")

# Loop over each source
for ((i=0; i<NUM_SOURCES; i++)); do
    # Read basic source info
    RUN_NAME=$(jq -r ".sources[$i].runName" "$CONFIG_FILE")
    SOURCE_NAME=$(jq -r ".sources[$i].source" "$CONFIG_FILE")
    DETECTOR=$(jq -r ".sources[$i].detector" "$CONFIG_FILE")

    # ---- Read C++ vectors for ROOT ----
    ENERGIES=$(jq -r ".sources[$i].energies | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    SIGMA_GUESS=$(jq -r ".sources[$i].sigma_guess | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    INTENSITY=$(jq -r ".sources[$i].intensity | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    DINTENSITY=$(jq -r ".sources[$i].dIntensity | map(tonumber) | join(\",\")" "$CONFIG_FILE")

    # ---- Read search windows ----
    # Convert JSON array into a Bash array
    SWIN_ARRAY=($(jq -r ".sources[$i].search_windows | map(tonumber) | join(\" \")" "$CONFIG_FILE"))

    # Construct pairs of search windows for ROOT
    SEARCH_WINDOWS="{"
    for ((j=0; j<${#SWIN_ARRAY[@]}; j+=2)); do
        w1=${SWIN_ARRAY[$j]}
        w2=${SWIN_ARRAY[$j+1]}
        SEARCH_WINDOWS+="{${w1},${w2}}"
        if [ $j -lt $((${#SWIN_ARRAY[@]}-2)) ]; then
            SEARCH_WINDOWS+=","
        fi
    done
    SEARCH_WINDOWS+="}"

    # ---- Read decay parameters ----
    HALFLIFE=$(jq -r ".sources[$i].halflife" "$CONFIG_FILE")
    D_HALFLIFE=$(jq -r ".sources[$i].dHalflife" "$CONFIG_FILE")
    ACTIVITY=$(jq -r ".sources[$i].activity" "$CONFIG_FILE")
    D_ACTIVITY=$(jq -r ".sources[$i].dActivity" "$CONFIG_FILE")
    INITIAL_ACTIVITY=$(jq -r ".sources[$i].initial_activity" "$CONFIG_FILE")
    MEASURE_TIME=$(jq -r ".sources[$i].measure_time" "$CONFIG_FILE")
    DECAY_TIME=$(jq -r ".sources[$i].decay_time" "$CONFIG_FILE")

    # ---- Display source info ----
    echo "----------------------------------------------------------"
    echo "Processing source $SOURCE_NAME"
    echo "  Run name: $RUN_NAME"
    echo "  Detector: $DETECTOR"
    echo "  Energies: $ENERGIES"
    echo "  Search windows: $SEARCH_WINDOWS"
    echo "  Sigma guess: $SIGMA_GUESS"
    echo "  Intensity: $INTENSITY ± $DINTENSITY"
    echo "  Half-life: $HALFLIFE ± $D_HALFLIFE"
    echo "  Activity: $ACTIVITY ± $D_ACTIVITY"
    echo "  Initial activity: $INITIAL_ACTIVITY"
    echo "  Measure time: $MEASURE_TIME"
    echo "  Decay time: $DECAY_TIME"
    echo "----------------------------------------------------------"

    # ----------------------------------------------------
    # 1. Run PeakIntegrator.C in ROOT
    # ----------------------------------------------------
    # Constructs a ROOT command to analyze peaks for this source
    ROOT_CMD=".L PeakIntegrator.C
peak_analysis(\"$RUN_NAME\", \"$DETECTOR\", \"$SOURCE_NAME\", std::vector<int>{$IDS}, SourcePeakData{
    \"$SOURCE_NAME\",
    std::vector<double>{$ENERGIES},
    $SEARCH_WINDOWS,
    std::vector<double>{$SIGMA_GUESS},
    std::vector<double>{$INTENSITY},
    std::vector<double>{$DINTENSITY},
    $HALFLIFE,
    $D_HALFLIFE,
    $ACTIVITY,
    $D_ACTIVITY,
    $INITIAL_ACTIVITY,
    $MEASURE_TIME,
    $DECAY_TIME
});"

    # Execute the ROOT command in batch mode
    echo "$ROOT_CMD"
    echo "$ROOT_CMD" | root -l -b

done

# ----------------------------------------------------
# 2. Combine output files if multiple sources exist
# ----------------------------------------------------
OUTPUT_FILES=()
for ((i=0; i<NUM_SOURCES; i++)); do
    DETECTOR=$(jq -r ".sources[$i].detector" "$CONFIG_FILE")
    SOURCE_NAME=$(jq -r ".sources[$i].source" "$CONFIG_FILE")
    OUTPUT_FILES+=("efficiency/PeakAreas_Detector=${DETECTOR}_Source=${SOURCE_NAME}.txt")
done

# ----------------------------------------------------
# 3. Run FitEfficiency.C
# ----------------------------------------------------
# If multiple sources exist, combine their peak area files into one
if [ ${#OUTPUT_FILES[@]} -gt 1 ]; then
    COMBINED_FILE="efficiency/PeakAreas_Detector=${DETECTOR}_Source=Combined.txt"
    # Take the header from the first file
    head -n 1 "${OUTPUT_FILES[0]}" > "$COMBINED_FILE"
    # Append all lines without headers from each file
    for f in "${OUTPUT_FILES[@]}"; do
        tail -n +2 "$f" >> "$COMBINED_FILE"
    done
    echo "✅ Files have been combined into $COMBINED_FILE"

    # Execute FitEfficiency.C on the combined file
    ROOT_CMD_2=".L FitEfficiency.C
FitEfficiency(\"$DETECTOR\", \"Combined\", $FIT_MODEL, $ENERGY_INTERPOLATION);"
    echo "$ROOT_CMD_2"
    echo "$ROOT_CMD_2" | root -l -b
fi

# ----------------------------------------------------
# Completion message
# ----------------------------------------------------
echo "=========================================================="
echo "✅ All efficiency calculations executed successfully."
echo "=========================================================="
