#!/bin/bash

# ====================================================================
# Cross Section SCRIPT
# This script automates the execution of CrossSection.C using a JSON
# configuration file containing all input parameters.
# ====================================================================

# Check if 'jq' is installed
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed. Please install it."
    echo "On Debian/Ubuntu, run: sudo apt-get install jq"
    exit 1
fi

# Check if JSON configuration file was provided
if [ -z "$1" ]; then
    echo "ERROR: Please provide the JSON configuration file as an argument."
    echo "Example: ./Run_CS.sh config_cs.json"
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

# =========================
# 1. Read parameters from JSON
# =========================

# Input Files
HIST_FILE=$(jq -r '.inputFiles.histogramFile' "$CONFIG_FILE")
HIST_NAME=$(jq -r '.inputFiles.histogramName' "$CONFIG_FILE")
GAMMA_ENERGY=$(jq -r '.inputFiles.gammaEnergy' "$CONFIG_FILE")

# Flux Info
FLUX_FILE_INPUT=$(jq -r '.fluxInfo.fluxFileInput' "$CONFIG_FILE")
FLUX_FILE_OUTPUT=$(jq -r '.fluxInfo.fluxFileOutput' "$CONFIG_FILE")
HIST_FLUX_NAME=$(jq -r '.fluxInfo.histogramFluxName' "$CONFIG_FILE")

# Rebinning & Pulses
REBIN_FACTOR=$(jq -r '.rebinningAndPulses.rebinFactor' "$CONFIG_FILE")
MIN_ENTRIES=$(jq -r '.rebinningAndPulses.minEntries' "$CONFIG_FILE")
NR_PULSES=$(jq -r '.rebinningAndPulses.nrPulses' "$CONFIG_FILE")

# Target Parameters
MASS_NUMBER=$(jq -r '.targetParameters.massNumber' "$CONFIG_FILE")
DENSITY=$(jq -r '.targetParameters.density' "$CONFIG_FILE")
DETECTOR_EFF=$(jq -r '.targetParameters.detectorEfficiency' "$CONFIG_FILE")

# Fit Parameters
PEAK_SIGMA=$(jq -r '.fitParameters.peakSearchSigma' "$CONFIG_FILE")
PEAK_RATIO=$(jq -r '.fitParameters.peakRatio' "$CONFIG_FILE")
PRE_FIT_RANGE=$(jq -r '.fitParameters.preFitRange' "$CONFIG_FILE")
GAMMA_WINDOW=$(jq -r '.fitParameters.gammaEnergyWindow' "$CONFIG_FILE")

# Output Files
PDF_NAME=$(jq -r '.outputFiles.pdfName' "$CONFIG_FILE")
ROOT_NAME=$(jq -r '.outputFiles.rootName' "$CONFIG_FILE")

# =========================
# 2. Display configuration
# =========================
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

echo "Executing ROOT macro..."
echo "$ROOT_CMD"
echo "$ROOT_CMD" | root -l -b

echo "=========================================================="
echo "✅ CrossSection.C execution completed."
echo "=========================================================="
