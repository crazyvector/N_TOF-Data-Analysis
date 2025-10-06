#!/bin/bash

# ====================================================================
# EFFICIENCY CALCULATION SCRIPT
# ====================================================================

# ----------------------------------------------------
# A. Verificări inițiale
# ----------------------------------------------------
# Check if 'jq' is installed
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

# Check if JSON configuration file was provided
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_Efficiency.sh config_eff.json"
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

# ----------------------------------------------------
# B. Parametri globali
# ----------------------------------------------------
FIT_MODEL=$(jq -r '.modelChoice' "$CONFIG_FILE")
ENERGY_INTERPOLATION=$(jq -r '.energyInterpolation' "$CONFIG_FILE")
IDS=$(jq -r '.ids | join(",")' "$CONFIG_FILE")

if [ "$FIT_MODEL" != "4" ] && [ "$FIT_MODEL" != "6" ]; then
    echo "⚠️ Modelul de fit ('modelChoice') nu este 4 sau 6. Se folosește 4."
    FIT_MODEL=4
fi

# ----------------------------------------------------
# C. Iterarea prin sursele de eficiență
# ----------------------------------------------------
NUM_SOURCES=$(jq -r '.sources | length' "$CONFIG_FILE")

for ((i=0; i<NUM_SOURCES; i++)); do
    RUN_NAME=$(jq -r ".sources[$i].runName" "$CONFIG_FILE")
    SOURCE_NAME=$(jq -r ".sources[$i].source" "$CONFIG_FILE")
    DETECTOR=$(jq -r ".sources[$i].detector" "$CONFIG_FILE")

    # Vectori C++
    ENERGIES=$(jq -r ".sources[$i].energies | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    SIGMA_GUESS=$(jq -r ".sources[$i].sigma_guess | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    INTENSITY=$(jq -r ".sources[$i].intensity | map(tonumber) | join(\",\")" "$CONFIG_FILE")
    DINTENSITY=$(jq -r ".sources[$i].dIntensity | map(tonumber) | join(\",\")" "$CONFIG_FILE")

    # Citim toate valorile ca un array
    SWIN_ARRAY=($(jq -r ".sources[$i].search_windows | map(tonumber) | join(\" \")" "$CONFIG_FILE"))

    # Construim perechile
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

    # Alte valori
    HALFLIFE=$(jq -r ".sources[$i].halflife" "$CONFIG_FILE")
    D_HALFLIFE=$(jq -r ".sources[$i].dHalflife" "$CONFIG_FILE")
    ACTIVITY=$(jq -r ".sources[$i].activity" "$CONFIG_FILE")
    D_ACTIVITY=$(jq -r ".sources[$i].dActivity" "$CONFIG_FILE")
    INITIAL_ACTIVITY=$(jq -r ".sources[$i].initial_activity" "$CONFIG_FILE")
    MEASURE_TIME=$(jq -r ".sources[$i].measure_time" "$CONFIG_FILE")
    DECAY_TIME=$(jq -r ".sources[$i].decay_time" "$CONFIG_FILE")

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
    # 1. Run Peak Integrator
    # ----------------------------------------------------
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

    echo "$ROOT_CMD"
    echo "$ROOT_CMD" | root -l -b

    # ----------------------------------------------------
    # 2. Run Efficiency Fit
    # ----------------------------------------------------
    ROOT_CMD_2=".L FitEfficiency.C
FitEfficiency(\"$DETECTOR\", \"$SOURCE_NAME\", $FIT_MODEL, $ENERGY_INTERPOLATION);"

    echo "$ROOT_CMD_2"
    echo "$ROOT_CMD_2" | root -l -b

done

echo "=========================================================="
echo "✅ All efficiency calculations executed successfully."
echo "=========================================================="
