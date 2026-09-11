# Changelog

All notable changes to this project should be documented here.


## 0.3.0 — 2026-09-11

- Added `rain_temp_timeseries.py` for time-series analysis of daily bias-adjusted precipitation and temperature data.

- Added monthly, seasonal, and annual aggregation of `prAdjust`, `tasmaxAdjust`, and `tasminAdjust`.

- Added derived mean temperature (`tasmeanAdjust`) calculated from daily maximum and minimum temperature.

- Added ensemble time-series uncertainty analysis using minimum, mean, and maximum across available GCM × RCM ensemble members.

- Added NetCDF outputs for gridded ensemble minimum, mean, and maximum time series.

- Added CSV outputs containing study-area ensemble minimum, mean, and maximum time series.

- Added rainfall and temperature uncertainty plots showing individual ensemble-member time series, ensemble min-max range, and ensemble mean.

- Added automatic discovery of downloaded scenarios, time windows, GCMs, RCMs, and study-area feature IDs from existing NARCliM workflow outputs.

- Added handling of model-specific climate calendars through common monthly, seasonal, and annual analysis periods.

- Added `flood.py` for approximate climate-adjusted flood-hazard analysis using Bureau of Meteorology IFD rainfall data and ARR climate-change factors.

- Added flood-hazard outputs for baseline and future climate scenarios and time horizons.

- Updated the package public API to expose `run_rain_temp_timeseries` and `run_flood_hazard`.


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