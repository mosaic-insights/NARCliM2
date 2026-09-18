NARCliM2 Workflow

A reusable Python package for discovering, downloading, spatially
subsetting, processing, summarising, and mapping NARCliM2.0 climate
data from the NCI THREDDS service. It supports climate-model ensemble
analysis, rainfall/temperature time series, storm indicators, spatial
summaries, and approximate climate-adjusted flood-frequency analysis.

Developer: Jabbar Khaledi --- Data and Geospatial Analyst | Python
Developer

Key capabilities

Access NARCliM2.0 through NCI THREDDS/OPeNDAP and automatically
select the appropriate spatial domain.

Extract climate data for polygon or point study areas.

Process multiple variables, scenarios, GCMs, RCMs, and user-defined
time horizons.

Calculate temporal-average rasters for individual GCM × RCM members.

Calculate pixel-wise ensemble minimum and maximum, with mean
optional.

Calculate climate change as future scenario/horizon minus historical
baseline.

Derive configurable annual storm-wind exceedance indicators from
sfcWindmax.

Generate monthly, seasonal, and annual rainfall/temperature time
series.

Summarise ensemble outputs for arbitrary polygon and point datasets.

Produce GeoTIFF, NetCDF, GeoPackage, CSV, plots, maps, and
processing manifests.

Run approximate climate-adjusted flood-frequency analysis using BoM
and ARR information.

Flood note: flood.py is a climate-risk screening/comparative
method, not a hydraulic flood model. It does not calculate flood
extent, depth, velocity, routing, drainage-network capacity,
floodplain storage, or property-scale inundation.

Installation

Clone the repository:

git clone https://github.com/mosaic-insights/NARCliM2.git
cd NARCliM2

Install the package in editable mode from the repository root (the
folder containing pyproject.toml):

python -m pip install -e .

Editable installation is recommended for development because changes
under src/narclim_workflow/ are available directly without
reinstalling the package.

For development dependencies:

python -m pip install -e ".[dev]"

Python 3.11+ is recommended. To confirm the installation:

import narclim_workflow
print(narclim_workflow.__file__)

NARCliM configuration

Domains

Domain                    Approx. resolution Coverage

NARCliM2-0-SEAus-04                 ~4 km South-East Australia
AUS-18                             ~18 km Australia

Use domain_key="auto" for automatic selection, or seaus_4km /
aus_18 to force a domain.

Climate models

GCMs: ACCESS-ESM1-5, EC-Earth3-Veg, MPI-ESM1-2-HR,
NorESM2-MM, UKESM1-0-LL

RCMs: NARCliM2-0-WRF412R3, NARCliM2-0-WRF412R5

With 5 GCMs × 2 RCMs, up to 10 ensemble members can contribute to
ensemble statistics.

Scenarios

historical, ssp126, ssp245, ssp370

Availability depends on variable, domain, scenario, GCM, RCM, and
period. Unavailable combinations are skipped and reported.

Climate variables

Variable                Description             Frequency

prAdjust              Bias-adjusted           daily
precipitation

tasmaxAdjust          Bias-adjusted maximum   daily
temperature

tasminAdjust          Bias-adjusted minimum   daily
temperature

TXge35                Days with maximum       yearly
temperature ≥35°C

TNlt2                 Days with minimum       yearly
temperature <2°C

FFDIgt50              Days with FFDI ≥50      yearly

R20mm                 Days with precipitation yearly
≥20 mm

R99p                  Total precipitation     yearly
from extremely wet days

SPI12                 12-month Standardised   monthly
Precipitation Index

tasmeanAdjust is derived as (tasmaxAdjust + tasminAdjust) / 2; it is
not downloaded separately.

StormDays is derived from daily sfcWindmax. For example,
storm_threshold_kmh=89.0 produces StormDaysGT89. The original daily
wind subset can optionally be retained for QA.

Workflow

Study area
    ↓
workflow.py
Download + spatially subset NARCliM data
    ↓
uncertainty.py
Temporal average for each GCM × RCM
    ↓
Pixel-wise ensemble Min / [optional Mean] / Max
    ↓
Scenario − historical baseline change
    ↓
spatial_summary.py
Summarise ensemble values and changes for polygons/points
    ↓
ensemble_maps.py
Create maps

For each variable × scenario × horizon, uncertainty.py first averages
each GCM × RCM member through time and then calculates ensemble
statistics across those member-average rasters.

change_min  = scenario_min  − baseline_min
change_max  = scenario_max  − baseline_max
change_mean = scenario_mean − baseline_mean   # when mean is enabled

For polygons, spatial_summary.py calculates the spatial mean of
intersecting cells from each ensemble-statistic raster. For points, it
samples the intersecting raster cell.

1. Download NARCliM data

from pathlib import Path
from narclim_workflow import run_workflow

VECTOR_PATH = Path(r"C:\path\to\study_area.shp")
OUTPUT_ROOT = Path(r"C:\NARCliM_Outputs")

TIME_WINDOWS = {
    "baseline": (1985, 2014),
    "short_term": (2021, 2040),
    "mid_term": (2041, 2060),
    "long_term": (2081, 2100),
}

manifest = run_workflow(
    vector_path=VECTOR_PATH,
    output_root=OUTPUT_ROOT,
    time_windows=TIME_WINDOWS,
    variables=["TXge35", "R20mm", "FFDIgt50", "SPI12", "StormDays"],
    scenarios=["historical", "ssp245", "ssp370"],
    gcms=[
        "ACCESS-ESM1-5", "EC-Earth3-Veg", "MPI-ESM1-2-HR",
        "NorESM2-MM", "UKESM1-0-LL",
    ],
    rcms=["NARCliM2-0-WRF412R3", "NARCliM2-0-WRF412R5"],
    domain_key="auto",
    id_field=None,
    overwrite=False,
    storm_threshold_kmh=89.0,
    save_storm_source=True,
)

print(manifest["status"].value_counts(dropna=False))

2. Ensemble uncertainty and climate change

from pathlib import Path
from narclim_workflow.uncertainty import run_uncertainty_workflow

INPUT_ROOT = Path(r"C:\NARCliM_Outputs")
BOUNDARY_PATH = Path(r"C:\path\to\study_area.shp")

manifest = run_uncertainty_workflow(
    input_root=INPUT_ROOT,
    boundary_path=BOUNDARY_PATH,
    variables=["TXge35", "R20mm", "FFDIgt50", "SPI12", "StormDaysGT89"],
    scenarios=["historical", "ssp126", "ssp245", "ssp370"],
    time_windows=["baseline", "short_term", "mid_term", "long_term"],
    gcms=[
        "ACCESS-ESM1-5", "EC-Earth3-Veg", "MPI-ESM1-2-HR",
        "NorESM2-MM", "UKESM1-0-LL",
    ],
    rcms=["NARCliM2-0-WRF412R3", "NARCliM2-0-WRF412R5"],
    feature_ids=None,
    boundary_id_field=None,
    output_folder_name="Climate_Indicies",
    include_mean=False,
    calculate_change=True,
    baseline_scenario="historical",
    baseline_time_horizon="baseline",
    overwrite=False,
    verbose=True,
)

include_mean=False produces Min/Max without Mean. Set it to True
when the ensemble mean is required.

3. Spatial summary

from pathlib import Path
from narclim_workflow.spatial_summary import run_spatial_summary

FEATURE_PATH = Path(r"C:\path\to\planning_zones.shp")
CLIMATE_INDICES_ROOT = Path(r"C:\NARCliM_Outputs\Climate_Indicies")

outputs = run_spatial_summary(
    feature_path=FEATURE_PATH,
    climate_indices_root=CLIMATE_INDICES_ROOT,
    dataset_name="Planning_Zones",
    feature_id_field=None,
    all_touched=True,
    include_mean=False,
    calculate_change=True,
    baseline_scenario="historical",
    baseline_time_horizon="baseline",
    overwrite=True,
    verbose=True,
)

Optional filters include variables, scenarios, time_horizons, and
source feature IDs.

4. Rainfall and temperature time series

rain_temp_timeseries.py preserves the time dimension rather than
averaging the whole horizon. Rainfall is aggregated as totals;
temperature as means. Supported frequencies are Monthly, Seasonal
(DJF/MAM/JJA/SON), and Annual.

from pathlib import Path
from narclim_workflow import run_rain_temp_timeseries

manifest = run_rain_temp_timeseries(
    input_root=Path(r"C:\NARCliM_Outputs"),
    variables=["prAdjust", "tasmaxAdjust", "tasminAdjust", "tasmeanAdjust"],
    frequencies=["Monthly", "Seasonal", "Annual"],
    scenarios=None,
    time_windows=None,
    output_folder_name="Climate_Indicies",
    overwrite=False,
    plot_dpi=300,
    show_plots=False,
    verbose=True,
)

Scenarios and time windows are automatically discovered when set to
None. Outputs include gridded NetCDF time series, CSV summaries, and
uncertainty plots.

5. Approximate flood-hazard analysis

flood.py uses BoM design rainfall and ARR climate-change/loss
information to estimate climate-driven changes in equivalent
rainfall-runoff event frequency. Locations can be coordinates or a point
shapefile/GeoPackage and can be identified as Urban or Rural.

from pathlib import Path
from narclim_workflow import run_flood_hazard

results = run_flood_hazard(
    output_root=Path(r"C:\NARCliM_Outputs"),
    time_windows={
        "baseline": (1985, 2014),
        "mid_term": (2041, 2060),
        "long_term": (2081, 2100),
    },
    scenarios=["historical", "ssp245", "ssp370"],
    location_path=Path(r"C:\path\to\flood_locations.gpkg"),
    name_field="Name",
    context_field="Context",
    overwrite=False,
    verbose=True,
)

Default reference ARIs are 1, 10, and 100 years. Default durations are 1
hour for Urban and 24 hours for Rural. Principal outputs are
Flood_Hazard.gpkg, Flood_Hazard.csv, and Flood_Runoff_Curves.csv.

6. Ensemble maps

from pathlib import Path
from narclim_workflow.ensemble_maps import plot_ensemble_maps

maps = plot_ensemble_maps(
    input_root=Path(r"C:\NARCliM_Outputs"),
    boundary_path=Path(r"C:\path\to\study_area.shp"),
    variables=None,
    scenarios=None,
    statistics=["min", "max"],  # add "mean" if generated
    cmap="viridis",
    dpi=300,
    overwrite=True,
    show=False,
    verbose=True,
)

Output structure

<OUTPUT_ROOT>/
├── NARCliM Data/
│   ├── <climate variables>/
│   ├── sfcWindmax/
│   ├── StormDaysGT<threshold>/
│   └── Flood/
└── Climate_Indicies/
    ├── Temporal_Average/
    ├── Ensemble_Stats/
    ├── Hazards/
    │   ├── Ensemble_summary_all_variables.gpkg
    │   ├── Ensemble_summary_<variable>.gpkg
    │   ├── <dataset>_Hazards.gpkg
    │   ├── <dataset>_Hazards.csv
    │   └── Flood/
    ├── Rain_Temp_Timeseries/
    └── Plot and Maps/

Climate_Indicies is retained because it is the current package
output-folder name. If this spelling is standardised later, update
package configuration, examples, and downstream paths together.

Interpretation and QA

Ensemble min/max: model-ensemble bounds after temporal
averaging, not minimum/maximum individual years within the horizon.

Polygon summaries: spatial means of the corresponding ensemble
rasters, not minimum/maximum cells inside the polygon.

Point summaries: values sampled from the intersecting raster
cell.

Time-series outputs: ensemble statistics are calculated at each
monthly, seasonal, or annual time step.

Raster QA: check CRS, extent, resolution, valid-cell count, and
grid alignment before project use.

Repository structure

NARCliM2/
├── README.md
├── pyproject.toml
├── example.py
├── NARCliM_Example.ipynb
├── CITATION.cff
├── CHANGELOG.md
├── .gitignore
└── src/
    └── narclim_workflow/
        ├── __init__.py
        ├── catalog.py
        ├── config.py
        ├── domain.py
        ├── spatial.py
        ├── workflow.py
        ├── uncertainty.py
        ├── spatial_summary.py
        ├── ensemble_maps.py
        ├── rain_temp_timeseries.py
        └── flood.py

For a public release, add an appropriate open-source LICENSE after
confirming the intended licensing terms.

Citation

If you use or adapt this package, please cite the package and developer.
A CITATION.cff file is included for GitHub citation metadata.

Jabbar Khaledi
Data and Geospatial Analyst | Python Developer
Email: jabbarkhaledi88@gmail.com
Primary language: Python