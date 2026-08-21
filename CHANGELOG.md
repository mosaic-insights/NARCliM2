# Changelog

All notable changes to this project should be documented here.

## 0.2.0 — 2026-08-21

- Added configurable `StormDays` processing from daily NARCliM `sfcWindmax`.
- Added user-defined `storm_threshold_kmh`.
- Added optional retention of original daily `sfcWindmax` study-area subsets.
- Added dynamic derived names such as `StormDaysGT89` and `StormDaysGT63`.
- Preserved three-decimal precision for rare `StormDaysGT*` frequencies.
- Standardised ensemble raster names, for example `TXge35_Mean_ssp245_2041_2060.tif`.
- Added spatial-summary CSV output alongside GeoPackage output.
- Added `time_horizon_year` to spatial-summary outputs.
- Added ensemble raster plotting with study-area boundary overlay.
- Expanded README and example documentation to cover the complete notebook workflow.
- Added GitHub support files (`.gitignore`, `CITATION.cff`) and development dependencies.

## 0.1.0

- Initial NARCliM2.0 THREDDS discovery, spatial extraction, and download workflow.
- Added climate-model uncertainty processing and ensemble statistics.
