# Efficiency Calculation Script

This script automates the calculation of detector efficiency for multiple sources, reading parameters from a JSON configuration file and executing ROOT macros (`PeakIntegrator.C` and `FitEfficiency.C`). It is designed to:

- Read global and source-specific parameters from a JSON file.
- Run a peak integration macro for each source.
- Combine output files if multiple sources exist.
- Fit efficiency curves using `FitEfficiency.C`.

---

## Workflow Description

1. **Initial Checks**
   - Verify that `jq` is installed. This tool is required to parse JSON configuration files. If `jq` is missing, the script exits with an error.
   - Verify that a JSON configuration file is provided as the first argument. If missing or if the file does not exist, the script exits.

2. **Read Global Parameters**
   - `FIT_MODEL`: The chosen fit model (usually 4 or 6). Defaults to 4 if an invalid model is provided.
   - `ENERGY_INTERPOLATION`: Defines the type of energy interpolation.
   - `IDS`: Detector IDs as a comma-separated string for use in ROOT macros.

3. **Iterate Through Sources**
   - For each source defined in the JSON file:
     - Read source information: run name, source name, and detector.
     - Read C++ vectors for ROOT: energies, sigma guesses, intensities, uncertainties.
     - Convert search windows from JSON array to ROOT-friendly pairs of `{min,max}`.
     - Read decay parameters: half-life, uncertainty in half-life, activity, uncertainty in activity, initial activity, measurement time, and decay time.
     - Display all information for verification.
     - Build the ROOT command to execute `PeakIntegrator.C`, passing all parameters for the current source.
     - Execute the ROOT macro using `root -l -b`.

4. **Combine Output Files**
   - If there are multiple sources, combine the `PeakAreas` output files into a single file:
     - Take the header from the first file.
     - Append all subsequent lines from all files, skipping headers.
   - This produces a combined file named like `PeakAreas_Detector=<DETECTOR>_Source=Combined.txt`.

5. **Run Efficiency Fit**
   - Use the combined file (or single-source file) as input for the `FitEfficiency.C` macro.
   - Execute the ROOT macro, passing the detector, combined source name, fit model, and energy interpolation.

6. **Final Output**
   - The script prints confirmation messages after each major step.
   - All ROOT macros are executed in batch mode (`-b`), and the console displays the commands being run.
   - At the end, the script prints a message confirming that all efficiency calculations have been successfully executed.

---

## Summary

This workflow ensures a fully automated process for calculating detector efficiency:

- Parameters are read from a JSON file.
- Peak integration is done for each source.
- Output files are combined if needed.
- Efficiency fitting is performed using the specified model and interpolation method.
- All steps are logged to the console for monitoring and debugging.
