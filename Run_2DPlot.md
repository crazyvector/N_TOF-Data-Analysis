# 2D Plot Bash Script Documentation

## Overview

This bash script, `Run_2DPlot.sh`, is designed to automate the creation of **2D histograms** (Neutron Energy vs Gamma Energy) using **ROOT macros**. It reads parameters from a **JSON configuration file** and runs the ROOT macro `Plot2D.C` in **batch mode**.

The script allows users to quickly generate histograms for multiple ROOT files, detectors, and run numbers without manually editing the macro each time.

---

## Dependencies

The script requires the following tools:

1. **Bash shell** (default on Linux/macOS)
2. **ROOT** (https://root.cern/)
   - The script runs ROOT in batch mode (`root -l -b`) and executes `Plot2D.C`.
3. **jq** (https://stedolan.github.io/jq/)
   - Command-line JSON processor to read the configuration file.

**Install `jq` on Debian/Ubuntu:**
```bash
sudo apt-get install jq
```

---

## Usage

### Running the Script
```bash
./Run_2DPlot.sh config_2D.json
```

- `config_2D.json`: JSON file containing all configuration parameters.
- The script will read this file, build the ROOT command, and execute the macro in batch mode.

---

## How It Works

### 1. Dependency Check

The script first checks if `jq` is installed.  
If not, it exits with an error and provides installation instructions.

```bash
if ! command -v jq &> /dev/null; then
    echo "ERROR: 'jq' is not installed."
    exit 1
fi
```

### 2. Input Validation

The script checks that a JSON file is provided as an argument and exists in the filesystem.  
If not, it exits with an error.

```bash
CONFIG_FILE="$1"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: JSON config file '$CONFIG_FILE' not found!"
    exit 1
fi
```

### 3. Reading JSON Configuration

The script extracts parameters from the JSON file using `jq`. Example:

```bash
ROOT_FOLDER=$(jq -r '.rootFile' "$CONFIG_FILE")
CALIB_FILE=$(jq -r '.calibrationFile' "$CONFIG_FILE")
DETECTOR=$(jq -r '.detectorName' "$CONFIG_FILE")
DETECTOR_COUNT=$(jq -r '.detectorCount' "$CONFIG_FILE")
```

**Extracted information includes:**

- **ROOT_FOLDER:** Path containing ROOT files for the runs.
- **CALIB_FILE:** Path to the calibration parameters file.
- **DETECTOR:** Detector tree name inside ROOT files.
- **DETECTOR_COUNT:** Total number of detectors.
- **SOURCE:** Source name.
- **FIT_TYPE:** Type of energy calibration fit.
- **PS_INT_THRESHOLD:** Minimum/Maximum pulse intensity to consider.
- **CONDITION:** Condition for pulse intensity to consider
- **DISTANCE:** Distance parameter for neutron energy calculation.
- **X/Y axis parameters:** Binning and range for histograms.
- **RUNS_ARRAY:** List of run numbers to include.

### 4. Constructing ROOT File Paths

The script iterates over the runs array and builds the full ROOT file paths:

```bash
for RUN in "${RUNS_ARRAY[@]}"; do
    ROOT_FILES+=("$ROOT_FOLDER/run$RUN.root")
done
```

### 5. Printing Configuration

All configuration values are printed for verification:

```bash
echo "ROOT folder: $ROOT_FOLDER"
echo "Calibration file: $CALIB_FILE"
echo "Detector: $DETECTOR ($DETECTOR_COUNT)"
```

This helps to ensure the script is reading the JSON correctly.

### 6. Constructing the ROOT Command

The script builds a ROOT command that:

1. Loads the macro `Plot2D.C` with compilation (`.L Plot2D.C+`)
2. Calls the `plot2D()` function with all parameters from JSON.
3. Passes the runs array in ROOT format `{ "run1", "run2", ... }`.

```bash
ROOT_CMD=".L Plot2D.C+
plot2D({$(printf '"%s",' "${RUNS_ARRAY[@]}" | sed 's/,$//')}, \"$DETECTOR\", $DETECTOR_COUNT, \"$SOURCE\", \"$FIT_TYPE\", ...);"
```

### 7. Executing the ROOT Command

The command is piped directly to ROOT in **batch mode**:

```bash
echo "$ROOT_CMD" | root -l -b
```

- `-l`: no splash screen
- `-b`: batch mode (no GUI)
- This executes the macro and generates 2D histograms as ROOT files and optionally PNG images.

---

## Output

After execution, the script produces:

1. **ROOT files** containing the histograms (`histograma_<Detector>_<Source>.root`)
2. **Optional PNG images** of 2D histograms stored in `2D_Plots/`
3. **Terminal output** showing progress, configuration, and any warnings.

---

## Summary

This bash script automates the generation of 2D histograms for multiple runs and detectors.  

- **Input:** JSON configuration file  
- **Process:** Reads config → Builds ROOT file list → Runs ROOT macro  
- **Output:** ROOT files + PNG images of histograms  

It allows physics users to efficiently analyze multiple datasets without modifying the ROOT macro manually.
