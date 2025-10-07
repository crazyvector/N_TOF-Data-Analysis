# Cross Section Workflow

This document explains the workflow of the **Cross Section Bash Script**, which automates the execution of `CrossSection.C` using a JSON configuration file.

---

## Script Overview

The script automates the calculation of nuclear cross sections using ROOT. It reads all necessary input parameters from a JSON file, sets up the execution environment, and runs the ROOT macro in batch mode.

---

## Step-by-Step Workflow

### 1. Dependency and Configuration Check
- Checks if `jq` is installed. `jq` is required for parsing JSON files.
- Ensures that a JSON configuration file is provided as a command-line argument.
- Verifies that the specified JSON file exists.

### 2. Read Parameters from JSON
Parameters are organized into sections in the JSON file:

#### Input Files
- `HIST_FILE`: Path to the input ROOT histogram file.
- `HIST_NAME`: Name of the histogram inside the ROOT file.
- `GAMMA_ENERGY`: Energy of the gamma line of interest.

#### Flux Information
- `FLUX_FILE_INPUT`: Input flux file (e.g., from a monitor detector).
- `FLUX_FILE_OUTPUT`: Output flux file after processing.
- `HIST_FLUX_NAME`: Name of the flux histogram.

#### Rebinning & Pulses
- `REBIN_FACTOR`: Factor to rebin histograms for better statistics.
- `MIN_ENTRIES`: Minimum number of entries required per bin.
- `NR_PULSES`: Number of pulses to consider in the analysis.

#### Target Parameters
- `MASS_NUMBER`: Mass number (A) of the target nucleus.
- `DENSITY`: Density of the target material.
- `DETECTOR_EFF`: Efficiency of the gamma detector.

#### Fit Parameters
- `PEAK_SIGMA`: Standard deviation used in peak search.
- `PEAK_RATIO`: Peak ratio parameter for fitting.
- `PRE_FIT_RANGE`: Range used for pre-fit of the peaks.
- `GAMMA_WINDOW`: Energy window for gamma selection.

#### Output Files
- `PDF_NAME`: Name of the output PDF file for plots.
- `ROOT_NAME`: Name of the output ROOT file containing processed results.

### 3. Display Configuration
- Prints all parameters to the terminal for verification.
- Helps ensure correct paths and values before ROOT execution.

### 4. Execute ROOT Macro
- Constructs a ROOT command string (`ROOT_CMD`) that:
  - Loads `CrossSection.C`.
  - Calls `RunCrossSection()` with all parameters read from the JSON.
- Pipes the command to ROOT in batch mode (`root -l -b`), so no GUI opens.
- Ensures reproducible and automated execution for multiple configurations.

### 5. Completion Message
- Prints a final message indicating that the execution of `CrossSection.C` is completed successfully.

---

## Summary

The **Cross Section Bash Script**:

- Automates nuclear cross-section calculations.
- Uses a JSON configuration file for flexible input management.
- Runs the ROOT macro in batch mode for efficiency.
- Handles input files, detector information, flux data, histogram parameters, fit settings, and output files.
- Provides logging and progress information for transparency.

This script is ideal for repetitive analysis of multiple histograms or experimental runs where cross-section extraction needs to be consistent and reproducible.
