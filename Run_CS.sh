#!/bin/bash

# ====================================================================
# Cross Section SCRIPT
# This script automates the execution of CrossSection.C using a JSON
# configuration file containing all input parameters.
# ====================================================================

# -------------------------
# 0. Pre-checks
# -------------------------
# Check if 'jq' is installed
# 'jq' is a command-line JSON processor used to parse the JSON configuration file.
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

# Check if the user provided a JSON configuration file as the first argument
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_CS.sh config_cs.json"
    exit 1
fi

CONFIG_FILE="$1"

# Verify that the specified JSON configuration file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi

echo "=========================================================="
echo "Reading configuration from $CONFIG_FILE"
echo "=========================================================="

# =========================
# 1. Read parameters from JSON
# =========================
# Using 'jq', extract the necessary parameters from the JSON configuration file.

# ---- Input Files ----
# Histogram file path, histogram name, and the gamma energy to be used
HIST_FILE=$(jq -r '.inputFiles.histogramFile' "$CONFIG_FILE")
HIST_NAME=$(jq -r '.inputFiles.histogramName' "$CONFIG_FILE")
GAMMA_ENERGY=$(jq -r '.inputFiles.gammaEnergy' "$CONFIG_FILE")

# ---- Flux Info ----
# Input and output flux files and the name of the flux histogram
FLUX_FILE_INPUT=$(jq -r '.fluxInfo.fluxFileInput' "$CONFIG_FILE")
FLUX_FILE_OUTPUT=$(jq -r '.fluxInfo.fluxFileOutput' "$CONFIG_FILE")
HIST_FLUX_NAME=$(jq -r '.fluxInfo.histogramFluxName' "$CONFIG_FILE")

# ---- Rebinning & Pulses ----
# Parameters controlling histogram rebinning and pulse selection
REBIN_FACTOR=$(jq -r '.rebinningAndPulses.rebinFactor' "$CONFIG_FILE")   # Factor to rebin the histogram
MIN_ENTRIES=$(jq -r '.rebinningAndPulses.minEntries' "$CONFIG_FILE")    # Minimum entries per bin
NR_PULSES=$(jq -r '.rebinningAndPulses.nrPulses' "$CONFIG_FILE")        # Number of pulses to consider

# ---- Target Parameters ----
# Physical properties of the target material
MASS_NUMBER=$(jq -r '.targetParameters.massNumber' "$CONFIG_FILE")      # Mass number A of target nucleus
DENSITY=$(jq -r '.targetParameters.density' "$CONFIG_FILE")             # Target density
DETECTOR_EFF=$(jq -r '.targetParameters.detectorEfficiency' "$CONFIG_FILE")  # Detector efficiency factor

# ---- Fit Parameters ----
# Parameters used for peak finding and fitting
PEAK_SIGMA=$(jq -r '.fitParameters.peakSearchSigma' "$CONFIG_FILE")     # Gaussian sigma for peak search
PEAK_RATIO=$(jq -r '.fitParameters.peakRatio' "$CONFIG_FILE")           # Peak ratio threshold
PRE_FIT_RANGE=$(jq -r '.fitParameters.preFitRange' "$CONFIG_FILE")      # Range to perform pre-fit
GAMMA_WINDOW=$(jq -r '.fitParameters.gammaEnergyWindow' "$CONFIG_FILE") # Energy window around gamma peaks

# ---- Output Files ----
# Names of the output PDF and ROOT files
PDF_NAME=$(jq -r '.outputFiles.pdfName' "$CONFIG_FILE")
ROOT_NAME=$(jq -r '.outputFiles.rootName' "$CONFIG_FILE")

# =========================
# 2. Display configuration
# =========================
# Print all read parameters to the console for verification
echo "Histogram File:    $HIST_FILE"
echo "Histogram Name:    $HIST_NAME"
echo "Gamma Energy:      $GAMMA_ENERGY"
echo "Input Flux File:   $FLUX_FILE_INPUT"
echo "Output Flux File:  $FLUX_FILE_OUTPUT"
echo "Histogram Flux:    $HIST_FLUX_NAME"
echo "Rebin Factor:      $REBIN_FACTOR"
echo "Min Entries:       $MIN_ENTRIES"
echo "Nr Pulses:         $NR_PULSES"
echo "Mass Number:       $MASS_NUMBER"
echo "Density:           $DENSITY"
echo "Detector Eff.:     $DETECTOR_EFF"
echo "Peak Sigma:        $PEAK_SIGMA"
echo "Peak Ratio:        $PEAK_RATIO"
echo "Pre-Fit Range:     $PRE_FIT_RANGE"
echo "Gamma Window:      $GAMMA_WINDOW"
echo "Output PDF:        $PDF_NAME"
echo "Output ROOT:       $ROOT_NAME"
echo "=========================================================="

# =========================
# 3. Execute ROOT macro
# =========================
# Construct a ROOT command that:
# 1. Loads the CrossSection.C macro
# 2. Calls RunCrossSection() with all the parameters read from JSON
# 3. Runs ROOT in batch mode (-b) without splash screen (-l)
ROOT_CMD=".L CrossSection.C
RunCrossSection(
\"$HIST_FILE\", \"$HIST_NAME\",
\"$FLUX_FILE_INPUT\", \"$FLUX_FILE_OUTPUT\",
\"$HIST_FLUX_NAME\", $GAMMA_ENERGY,
$REBIN_FACTOR, $MIN_ENTRIES, $NR_PULSES,
$MASS_NUMBER, $DENSITY, $DETECTOR_EFF,
$PEAK_SIGMA, $PEAK_RATIO, $PRE_FIT_RANGE, $GAMMA_WINDOW,
\"$PDF_NAME\", \"$ROOT_NAME\");
"

# Print the constructed ROOT command to the console for debugging purposes
echo "Executing ROOT macro..."
echo "$ROOT_CMD"

# Pipe the ROOT command into ROOT itself in batch mode
echo "$ROOT_CMD" | root -l -b

# Completion message to indicate the script finished successfully
echo "=========================================================="
echo "✅ CrossSection.C execution completed."
echo "=========================================================="
