# ICESat-2 MIZ Width Analysis

This project provides a MATLAB-based analysis pipeline for calculating the Marginal Ice Zone (MIZ) width from ICESat-2 data. The code processes along-track sea ice data to identify the MIZ and compute its width based on sea ice concentration.

## Project Structure

The repository is organized into the following directories:

- **`Data/`**: Contains input data used in the analysis.
- **`Figure-IS2-Comp/`**: Scripts for generating figures from ICESat-2 data.
- **`Figure-SAR-comp/`**: Scripts for generating figures from SAR data.
- **`Geographic_Comp/`**: Scripts for geographic comparisons.
- **`compute_MIZ_width/`**: Core scripts for computing the MIZ width.
- **`Utilities/`**: Helper functions used across the analysis pipeline.

## How to Run the Analysis

The main driver script for the MIZ width computation is `compute_MIZ_width/drive_MIZ_width.m`. To run the analysis, follow these steps:

1. **Set up the environment**:
   - Ensure you have MATLAB installed.
   - Update the paths in `compute_MIZ_width/drive_MIZ_width.m` to point to your data and output directories. Specifically, modify the `OS_string`, `OPTS.track_folder`, and `OPTS.output_folder` variables.

2. **Run the main script**:
   - Open MATLAB and navigate to the project's root directory.
   - Run the main driver script:
     ```matlab
     run('compute_MIZ_width/drive_MIZ_width.m')
     ```

This will execute the full analysis pipeline, which includes:
- Preprocessing individual ICESat-2 tracks.
- Calculating along-track statistics.
- Referencing the statistics to the MIZ to determine its width.

The output files will be saved in the directory specified by `OPTS.output_folder`.