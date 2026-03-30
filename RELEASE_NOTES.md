# Release Notes

## v2.0 — 2026-03-30

This release represents a major reorganization and expansion of the codebase following the v1 archive. It is the version associated with the submitted/published manuscript.

### Repository restructuring

- Renamed `compute_MIZ_width/` → `Processing/` and consolidated all processing utilities under `Utilities/Processing/`
- Replaced the old `Figure-IS2-Comp/`, `Figure-SAR-comp/`, and `Geographic_Comp/` directories with a unified `Figures/` tree containing `Figure_1_PM_Comp/`, `Figures_IS2/`, and `Figure_SAR_comp/`
- Added `Drive_Figures.m` as a single master entry point for all figure generation
- Moved and reorganised figure utilities into `Utilities/Figures/`
- Removed `Geographic_Comp/compare_geographic_MIZ.m` (geographic comparison folded into main figure pipeline)
- Removed `Utilities/make_dummy_data.m`

### Processing pipeline

- Updated processing to jointly handle ICESat-2 data versions v6 and v7 (`AT_stats_SH_v6+7_all`)
- Added `nodark` photon classification (Linear Ice Fraction without dark returns: `LIF_nodark`)
- Improved performance of `generate_AT_statistics.m`; rewrote as a standalone utility under `Utilities/Analysis/`
- Added `reference_stats_to_MIZ.m` to `Utilities/Processing/` (previously at repo root)
- Removed the non-OSCAR OSCAR-SLA dependency; merged into the main branch

### Passive microwave comparison (Figure 1)

- Added full PM comparison pipeline: `calc_PM_stats.m`, `plot_PM_global_stats.m`, `plot_PM_bias_map.m`
- Coastal mask updated to **50 km** (was 25 km in v1)
- Added `create_coast_dist_data.m` for generating the coastal distance mask

### ICESat-2 figures

- Added `create_composite_figure.m` — main IS2 vs. PM bias comparison figure
- Added `create_parametric_plots.m` and supporting helpers (`get_param_plot_vectors.m`, `param_plot.m`)
- Added `fit_AMSR_offset.m` for inter-algorithm bias correction
- Rewrote `create_MIZ_width_panel.m` (expanded from 67 to 81 lines)
- Expanded `load_MIZ_waves.m` significantly (added wave classification logic, ~100 additional lines)
- Added `create_location_data.m`
- Updated track lists in `Track_Lists/` (trimmed `out_all.txt`)

### Supporting figures

- Added `Figures/Supporting-Figures/bias_plot_BS.m`
- Added presentation-quality figures under `Figures-For-Presentations/` including `make_global_figures_pres.m`, `create_composite_figure_pres.m`, `compare_phase_space.m`, and others

### Bug fixes and minor changes

- Fixed a bug in the MIZ edge detection (commit `922b144`)
- Small plot updates and refinements across multiple figure scripts

---

## v1.0 — 2025-03-20

Initial archive release. Core pipeline for computing along-track ICESat-2 MIZ statistics and referencing measurements to the MIZ edge. Basic figure scripts for IS2 and SAR comparisons.
