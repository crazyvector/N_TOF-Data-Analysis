# Energy Calibration Workflow

This document explains the workflow of the **Energy Calibration Bash Script**, which automates calibration for multiple runs, detectors, and sources using ROOT and JSON configuration files.

---

## Script Overview

The script handles four main calibration steps:

1. Summing or creating individual histograms per run.
2. Fitting gamma peaks for energy calibration.
3. Applying the selected fit type for all detectors.
4. Converting channel numbers to calibrated energy for each source.

It uses a JSON configuration file to define:

- Input/output folders
- Runs and sources
- Detector information
- Histogram and fit parameters

The workflow executes ROOT scripts in batch mode (`root -l -b`) for automation.

---

## Step-by-Step Workflow

### 1. Dependency and Configuration Check
- Checks if `jq` is installed for parsing JSON.
- Verifies that a configuration JSON file is provided and exists.

### 2. Global Parameter Reading
- Reads global parameters from the JSON file, including:
  - Input and output folders
  - Fit type
  - Histogram settings (rebin factor, x-axis min/max, bin width)
- Ensures folder paths end with a trailing slash.

### 3. Iteration Over Sources
For each source defined in the JSON:

- Extracts source-specific parameters:
  - Run names
  - Detector and number of detectors
  - Source name
  - IDs for fitting
  - Known gamma energies
  - Fit parameters (ratio, sigma, plot range, constants)
- Converts space-separated strings into Bash arrays for ROOT usage.

### 4. Calibration 1 – Histogram Creation
- Creates summed histograms for multiple runs or individual histograms for single runs.
- Uses the ROOT script `1_CalibrationsD.C`.
- Skips histogram creation if the output file already exists.

### 5. Calibration 2 – Peak Fitting
- Fits gamma peaks to extract energy calibration.
- Uses the ROOT script `2_CalibrationsD.C`.
- Passes detector IDs, energies, and fit parameters to ROOT.

### 6. Calibration 3 – Fit Type Selection
- Applies the chosen fit type (e.g., edge or Compton) to all detectors.
- Uses the ROOT script `3_CalibrationsD.C`.
- Writes results into an output file for later use.

### 7. Calibration 4 – Channel to Energy Conversion
- Converts channel numbers to calibrated energy for all runs and sources.
- Uses the ROOT script `4_CalibrationsD.C`.
- Iterates over sources again to perform the conversion.

### 8. Automation & Execution
- All ROOT scripts are executed in batch mode (`root -l -b`) without opening the GUI.
- The script prints progress messages and executed commands for transparency.

---

## Summary

The **Energy Calibration Bash Script** automates the entire calibration workflow for multiple runs, sources, and detectors. Its main advantages are:

- Fully automated ROOT execution.
- Flexible JSON-based configuration.
- Handles multiple sources and detectors.
- Skips already-processed histograms to save time.
- Maintains detailed logging of runs and parameters.

This script is ideal for repetitive calibration tasks in nuclear or particle physics experiments where multiple detectors and runs are involved.
