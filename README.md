# compare_MIZ_width

A MATLAB codebase for analyzing the width of the Southern Hemisphere **Marginal Ice Zone (MIZ)** using ICESat-2 (IS2) satellite lidar data, and for comparing sea ice concentration (SIC) estimates from multiple passive microwave (PM) satellite products. This repository supports the analysis and figures in an associated research manuscript.

---

## Overview

The MIZ is the transition zone between consolidated sea ice and open ocean. This project:

- Processes along-track ICESat-2 ATL07 sea ice segment data to characterize the MIZ
- Computes along-track statistics including Linear Ice Fraction (LIF), Wave-Affected Fraction (WAF), floe size distribution (FSD), and sea surface height (SSH)
- References all measurements to the MIZ edge (defined by the 80% sea ice concentration contour)
- Compares ICESat-2 MIZ width estimates against gridded passive microwave SIC products (AMSR2-NT, AMSR2-Bootstrap, NSIDC-CDR/SSMI)
- Validates results against co-located Sentinel-1 SAR classifications
- Produces publication-quality figures

---

## Repository Structure

```
compare_MIZ_width/
├── Processing/
│   ├── drive_MIZ_width.m           # Main driver: runs the full processing pipeline
│   ├── ANALYSIS_par.sh             # SLURM batch script for HPC execution
│   └── Utilities/
│       └── Processing/
│           ├── preprocess_single_track.m   # Reads and cleans individual IS2 HDF5 tracks
│           ├── create_AT_stat_file.m       # Loops over all tracks; parallel beam processing
│           ├── generate_AT_statistics.m    # Computes LIF, WAF, SSH, FSD per beam
│           └── reference_stats_to_MIZ.m   # Adds distance-to-MIZ-edge to each measurement
│
├── Figures/
│   ├── Drive_Figures.m             # Master figure driver; calls all sub-scripts in order
│   ├── Figure_1_PM_Comp/
│   │   ├── calc_PM_stats.m         # Loads and preprocesses gridded SIC datasets
│   │   ├── plot_PM_global_stats.m  # Time series and histograms of SIC product differences
│   │   └── plot_PM_bias_map.m      # Geographic bias maps between PM products
│   ├── Figures_IS2/
│   │   ├── load_MIZ_waves.m        # Loads MIZ data; classifies tracks by wave activity
│   │   ├── create_location_data.m  # Geographic context for track locations
│   │   ├── create_composite_figure.m   # Main IS2 vs. PM comparison figure
│   │   ├── create_MIZ_width_panel.m    # MIZ width analysis panel
│   │   ├── create_parametric_plots.m   # Parametric analysis figures
│   │   ├── fit_AMSR_offset.m           # Fits bias correction between AMSR products
│   │   └── Track_Lists/
│   │       ├── out_all.txt         # All processed tracks
│   │       ├── out_waves.txt       # Tracks with wave activity
│   │       ├── out_nowaves.txt     # Tracks without wave activity
│   │       ├── out_somewaves.txt   # Partially wave-affected tracks
│   │       └── out_wavy.txt        # Highly wave-affected subset
│   └── Figure_SAR_comp/
│       ├── comp_sentinel_IS2.m     # Co-located Sentinel-1 / IS2 comparison
│       └── make_sentinel_figure.m  # Generates SAR comparison figure
│
└── Data/                           # Reference data (coastal mask, etc.)
```

---

## Input Data

| Dataset | Description | Path (default) |
|---|---|---|
| ICESat-2 ATL07 | HDF5 sea ice segment files (v6+7) | `/gpfs/data/epscor/chorvat/IS2/Data/All_Track_Data/v6+7/SH/` |
| AMSR2-NT | AMSR2 NASA Team SIC grids | `/Data/SIC_Data/AMSR2-NT/` |
| AMSR2-Bootstrap | AMSR2 Bootstrap SIC grids | `/Data/SIC_Data/AMSR2-Bootstrap/` |
| NSIDC-CDR | CDR daily SIC (SSMI/SSMIS) | `/Data/SIC_Data/NSIDC-CDR/` |
| Sentinel-1 SAR | Classified SAR images with IS2 overlap list | — |
| Coastal mask | Distance-to-coast grid (50 km cutoff) | Bundled reference file |

All paths are configurable via `OPTS` in the driver scripts.

---

## Output

| File | Contents |
|---|---|
| `AT_stats_SH_v6+7_all.mat` | Full along-track statistics (`IS2_DATA`) and MIZ-referenced data (`MIZ_DATA`) for all beams and tracks |
| Figures | Publication figures saved to the configured `plot_save_str` directory |

### `MIZ_DATA` fields

| Field | Description |
|---|---|
| `lat`, `lon` | Geographic coordinates |
| `D_to_MIZ` | Distance to MIZ edge (m) |
| `H` | Height relative to sea surface |
| `E` | Height variance |
| `LIF` | Linear Ice Fraction |
| `WAF` | Wave-Affected Fraction |
| `SIC`, `SIC_amsr` | Sea ice concentration (CDR and AMSR2) |
| `RFSD`, `MFSD` | Floe size distribution metrics |

---

## Running the Code

### Local / interactive MATLAB

```matlab
% Step 1 — process all ICESat-2 tracks
run('Processing/drive_MIZ_width.m')

% Step 2 — generate all figures
run('Figures/Drive_Figures.m')
```

### HPC (SLURM)

```bash
sbatch Processing/ANALYSIS_par.sh
```

The batch script requests 10 hours, 16 cores, and 120 GB memory.

### Execution order

1. `drive_MIZ_width.m` — sets `OPTS`, calls `create_AT_stat_file` then `reference_stats_to_MIZ`
2. `Drive_Figures.m` — loads `.mat` output, then calls figure sub-scripts in sequence

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `OPTS.AT_window` | `[6250 6250]` m | Moving-window size for along-track statistics |
| `OPTS.AT_resolution` | `6250` m | Output bin resolution |
| `OPTS.do_weak` | `1` | Include weak laser beams (all 6 beams) |
| `SIC_threshold` | `0.80` | SIC contour used to define the MIZ edge |
| `wave_thresh` | `0.075` | WAF threshold for classifying wave-affected tracks |
| `cutoff_N` | `100` | Minimum valid measurements required per bin |
| Coastal mask | `50 km` | Minimum distance from coast for valid SIC data |

---

## Dependencies

- **MATLAB** (tested with R2021a+)
- MATLAB Mapping Toolbox
- MATLAB Statistics and Machine Learning Toolbox
- Plot-Tools — custom plotting utilities (`jbfill`, etc.)
- NE_Coastlines — coastline data for map figures

---

## Citation

If you use this code, please cite the associated manuscript (details to be added upon publication).

---

## License

This repository is released for archival and reproducibility purposes. Contact the authors for reuse permissions.
