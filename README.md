NARCliM2 Data and Hazards Analysis Workflow

A reusable Python workflow for NARCliM2.0 climate-data access,
processing, ensemble analysis, spatial summarisation, time-series
analysis, mapping, storm indicators, and climate-adjusted
flood-frequency screening.

The package is designed for reproducible climate-risk and geospatial
analysis across user-defined study areas, climate variables, model
ensembles, scenarios, and planning horizons.

Developer: Jabbar Khaledi --- Data and Geospatial Analyst |
Python Developer

Contents

Overview

Installation

NARCliM configuration

Package workflows

Core processing concepts

Examples

Output structure

QA and interpretation

Repository structure

Citation

Overview

Main capabilities

The package supports:

automated NARCliM2.0 access through NCI THREDDS/OPeNDAP and
catalogue discovery;

automatic selection of the appropriate NARCliM spatial domain;

spatial extraction for polygon and point study areas;

processing of multiple variables, scenarios, GCMs, RCMs, and time
horizons;

temporal averaging for individual GCM × RCM ensemble members;

pixel-wise ensemble Min / optional Mean / Max calculations;

future-climate change relative to historical baseline;

configurable storm-wind threshold indicators derived from
sfcWindmax;

retention of original daily storm-source data for QA when requested;

monthly, seasonal, and annual rainfall/temperature time-series
analysis;

derived daily mean temperature, tasmeanAdjust;

spatial summaries for planning zones, land use, catchments,
administrative areas, sites, assets, and other user datasets;

GeoTIFF, NetCDF, GeoPackage, CSV, plot, map, and processing-manifest
outputs;

approximate climate-adjusted flood-frequency analysis using BoM
design rainfall and ARR climate-change/loss information.

Design principles

Reproducibility --- climate-data processing is implemented as
repeatable Python workflows rather than manual GIS steps.
Flexibility --- study areas, variables, scenarios, time horizons,
models, statistics, and selected hazard parameters are configurable.
Transparency --- intermediate products, temporal averages, ensemble
statistics, metadata, source data, and manifests are retained where
useful.
Spatial applicability --- climate-model outputs are converted into
analysis-ready products for GIS, planning, and climate-risk assessment.

Installation

1. Clone the repository

git clone https://github.com/mosaic-insights/NARCliM2.git
cd NARCliM2

2. Install the package

Run the following from the repository root --- the folder containing
pyproject.toml:

python -m pip install -e .

The -e option installs the package in editable mode, which is
recommended for development because changes made under
src/narclim_workflow/ are available without reinstalling the package.

For development dependencies:

python -m pip install -e ".[dev]"

Python 3.11+ is recommended.

Optional installation check:

import narclim_workflow
print(narclim_workflow.__file__)

NARCliM configuration

Spatial domains

Domain                    Approx. resolution Coverage

NARCliM2-0-SEAus-04                 ~4 km South-East Australia
AUS-18                             ~18 km Australia

Domain selection is controlled with:

domain_key="auto"       # automatically select the appropriate domain
domain_key="seaus_4km"  # force ~4 km South-East Australia domain
domain_key="aus_18"     # force ~18 km Australia-wide domain

With domain_key="auto", the ~4 km domain is selected when the
supplied study area is fully within its stored footprint; otherwise the
workflow uses the Australia-wide domain.

Climate variables

Configured NARCliM source variables and climate indices are:

variables = [
    "prAdjust",       # daily bias-adjusted precipitation
    "tasmaxAdjust",   # daily bias-adjusted maximum temperature
    "tasminAdjust",   # daily bias-adjusted minimum temperature
    "TXge35",         # yearly number of days Tmax >= 35 C
    "TNlt2",          # yearly number of days Tmin < 2 C
    "FFDIgt50",       # yearly number of days with FFDI >= 50
    "R20mm",          # yearly number of days precipitation >= 20 mm
    "R99p",           # yearly precipitation from extremely wet days
    "SPI12",          # 12-month Standardised Precipitation Index
    "StormDays",      # derived annual wind-threshold exceedance indicator
]

Variable          Description       NCI branch               Frequency

prAdjust        Daily             bias-adjusted-output   day
bias-adjusted
precipitation

tasmaxAdjust    Daily             bias-adjusted-output   day
bias-adjusted
maximum
temperature

tasminAdjust    Daily             bias-adjusted-output   day
bias-adjusted
minimum
temperature

TXge35          Days with maximum bias-adjusted-output   yr
temperature ≥35°C

TNlt2           Days with minimum bias-adjusted-output   yr
temperature <2°C

FFDIgt50        Days with FFDI    DD                     yr
≥50

R20mm           Days with         DD                     yr
precipitation ≥20
mm

R99p            Total             DD                     yr
precipitation
from extremely
wet days

SPI12           12-month          DD                     mon
Standardised
Precipitation
Index

tasmeanAdjust is derived within rain_temp_timeseries.py:

tasmeanAdjust = (tasmaxAdjust + tasminAdjust) / 2

It is therefore not downloaded as a separate NARCliM source variable.

Storm processing

StormDays is also a derived indicator. The package reads daily
maximum near-surface wind speed, sfcWindmax, and applies a
user-defined threshold:

storm_threshold_kmh=89.0
save_storm_source=True

For example, 89 km/h produces:

StormDaysGT89

The indicator represents the annual number of days where daily
sfcWindmax exceeds the threshold. The original daily sfcWindmax
subset can be retained for QA and testing alternative thresholds without
downloading the source again.

Because high-wind exceedances can be rare, StormDaysGT* outputs retain
three decimal places; most other exported climate values use one
decimal place.

Scenarios

Configured scenarios are:

scenarios = [
    "historical",
    "ssp126",
    "ssp245",
    "ssp370",
]

Variables unavailable for a requested scenario/model combination are
skipped and reported rather than terminating the full batch. For
example, the package configuration treats SPI12 as unavailable for
historical.

Time horizons

Time windows are fully user configurable. A typical configuration is:

TIME_WINDOWS = {
    "baseline":   (1985, 2014),
    "short_term": (2021, 2040),
    "mid_term":   (2041, 2060),
    "long_term":  (2081, 2100),
}

Historical NARCliM data are generally within 1951–2014, while
projection scenarios cover 2015–2100. Requested windows should always
be checked against the datasets used in the analysis.

Climate models

Five GCMs are configured:

gcms = [
    "ACCESS-ESM1-5",
    "EC-Earth3-Veg",
    "MPI-ESM1-2-HR",
    "NorESM2-MM",
    "UKESM1-0-LL",
]

Two RCMs are configured:

rcms = [
    "NARCliM2-0-WRF412R3",
    "NARCliM2-0-WRF412R5",
]

With 5 GCMs × 2 RCMs, each grid cell can contain up to 10
ensemble-member values. Actual availability depends on variable,
domain, scenario, model, and period.

Package workflows

Core NARCliM workflow

Select study area
        ↓
Select domain automatically or explicitly
        ↓
workflow.py
        ↓
Discover NCI THREDDS catalogues
        ↓
Open remote NetCDF through OPeNDAP
        ↓
Subset requested time period and spatial extent
        ↓
Save local NARCliM NetCDF subsets
        ↓
uncertainty.py
        ↓
Temporal average for each GCM × RCM member
        ↓
Pixel-wise ensemble Min / [optional Mean] / Max
        ↓
Calculate scenario change relative to historical baseline
        ↓
spatial_summary.py
        ↓
Summarise climate values and changes for polygons or points
        ↓
ensemble_maps.py
        ↓
Create maps from ensemble rasters

Rainfall / temperature time-series workflow

This pathway preserves temporal variability instead of averaging the
complete planning horizon:

Daily prAdjust / tasmaxAdjust / tasminAdjust
        ↓
rain_temp_timeseries.py
        ↓
Derive tasmeanAdjust if requested
        ↓
Monthly / Seasonal / Annual aggregation
        ↓
Align ensemble-member time coordinates
        ↓
Pixel-wise ensemble Min / Mean / Max through time
        ↓
NetCDF + study-area CSV + uncertainty plots

Aggregation rules:

Rainfall:     Monthly total | Seasonal total | Annual total
Temperature:  Monthly mean  | Seasonal mean  | Annual mean
Seasons:      DJF | MAM | JJA | SON

December is assigned to the following DJF season year, incomplete
seasons are removed, and annual values are retained only where all 12
months are represented.

Flood-hazard workflow

Point locations / coordinates
        ↓
Identify Urban / Rural context
        ↓
BoM design rainfall + ARR climate/loss information
        ↓
flood.py
        ↓
Current rainfall-runoff reference events
        ↓
Climate-adjust future rainfall and relevant losses
        ↓
Future runoff-frequency curve
        ↓
Equivalent future flood frequency
        ↓
GeoPackage + CSV + runoff curves + QA/source files

flood.py is an approximate climate-driven rainfall-runoff frequency
method. It is not a hydraulic flood model and does not calculate
flood extent, depth, velocity, routing, drainage-network capacity,
floodplain storage, or property-scale inundation.

Core processing concepts

Temporal averaging

For each:

variable × scenario × GCM × RCM × grid cell × time horizon

the climate values are averaged across the selected horizon. The result
is one temporal-average raster for each GCM × RCM ensemble member.

Ensemble uncertainty

The temporal-average member rasters are combined and ensemble statistics
are calculated independently at each grid cell:

Min = lowest GCM × RCM temporal-average value
Mean = average GCM × RCM temporal-average value   [optional]
Max = highest GCM × RCM temporal-average value

Mean can be disabled:

include_mean=False

The Min/Max outputs are model-ensemble bounds, not the minimum and
maximum individual years within the horizon.

Climate change relative to baseline

When enabled:

calculate_change=True
baseline_scenario="historical"
baseline_time_horizon="baseline"

the package calculates:

change_min  = scenario_min  - baseline_min
change_max  = scenario_max  - baseline_max
change_mean = scenario_mean - baseline_mean   # only when Mean is enabled

This convention is used consistently by the uncertainty and
spatial-summary workflows.

Spatial summarisation

The ensemble products can be applied to planning zones, land use,
districts, LGAs, catchments, subcatchments, properties, sites, assets,
and other point/polygon datasets.

Polygon → spatial mean of intersecting cells from each ensemble raster
Point   → value sampled from the raster cell containing the point

Therefore, polygon fields named min and max are spatial means of the
ensemble-Min and ensemble-Max rasters; they are not the minimum and
maximum raster cells inside the polygon.

Examples

1. Download and subset NARCliM data

from pathlib import Path
from narclim_workflow import run_workflow

VECTOR_PATH = Path(r"C:\path\to\study_area.shp")
OUTPUT_ROOT = Path(r"C:\NARCliM_Outputs")

TIME_WINDOWS = {
    "baseline":   (1985, 2014),
    "short_term": (2021, 2040),
    "mid_term":   (2041, 2060),
    "long_term":  (2081, 2100),
}

manifest = run_workflow(
    vector_path=VECTOR_PATH,
    output_root=OUTPUT_ROOT,
    time_windows=TIME_WINDOWS,
    variables=[
        "prAdjust", "tasmaxAdjust", "tasminAdjust",
        "TXge35", "R99p", "FFDIgt50", "SPI12", "StormDays",
    ],
    scenarios=["historical", "ssp245", "ssp370"],
    gcms=[
        "ACCESS-ESM1-5", "EC-Earth3-Veg", "MPI-ESM1-2-HR",
        "NorESM2-MM", "UKESM1-0-LL",
    ],
    rcms=["NARCliM2-0-WRF412R3", "NARCliM2-0-WRF412R5"],
    domain_key="auto",
    id_field=None,
    overwrite=False,
    print_traceback=True,
    storm_threshold_kmh=89.0,
    save_storm_source=True,
)

print(manifest["status"].value_counts(dropna=False))

Large daily-variable runs can be substantial. Keep overwrite=False
when existing outputs should be reused.

2. Climate-model uncertainty and change

For a storm threshold of 89 km/h, use the derived variable name
StormDaysGT89:

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

    include_mean=False,                 # True to also calculate ensemble Mean
    calculate_change=True,
    baseline_scenario="historical",
    baseline_time_horizon="baseline",

    overwrite=False,
    verbose=True,
)

print(manifest["status"].value_counts(dropna=False))

Typical ensemble outputs include:

TXge35_Min_ssp245_2041_2060.tif
TXge35_Max_ssp245_2041_2060.tif
TXge35_Change_Min_ssp245_2041_2060.tif
TXge35_Change_Max_ssp245_2041_2060.tif

When include_mean=True, Mean and Change-Mean products are also
generated.

3. Spatial summary for planning zones or other features

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

    include_mean=False,                 # True if Mean rasters were generated
    calculate_change=True,
    baseline_scenario="historical",
    baseline_time_horizon="baseline",

    overwrite=True,
    verbose=True,
)

To restrict processing:

outputs = run_spatial_summary(
    feature_path=FEATURE_PATH,
    climate_indices_root=CLIMATE_INDICES_ROOT,
    variables=["TXge35", "FFDIgt50"],
    scenarios=["ssp245", "ssp370"],
    time_horizons=["mid_term", "long_term"],
    include_mean=False,
    calculate_change=True,
    all_touched=True,
    overwrite=True,
)

Typical outputs:

Planning_Zones_Hazards.gpkg
Planning_Zones_Hazards.csv

The GeoPackage contains climate-variable layers and the feature-based
change summary produced by the current spatial-summary workflow.

4. Ensemble maps

from pathlib import Path
from narclim_workflow.ensemble_maps import plot_ensemble_maps

maps = plot_ensemble_maps(
    input_root=Path(r"C:\NARCliM_Outputs"),
    boundary_path=Path(r"C:\path\to\study_area.shp"),
    variables=None,
    scenarios=None,
    statistics=["min", "max"],  # add "mean" when Mean rasters exist
    cmap="viridis",
    dpi=300,
    overwrite=True,
    show=False,
    verbose=True,
)

print(f"Created/found {len(maps)} map(s).")

Maps are written under:

Climate_Indicies/Plot and Maps/

5. Rainfall and temperature time-series analysis

Daily prAdjust, tasmaxAdjust, and tasminAdjust subsets must
already exist under NARCliM Data, normally from run_workflow().

from pathlib import Path
from narclim_workflow import run_rain_temp_timeseries

manifest = run_rain_temp_timeseries(
    input_root=Path(r"C:\NARCliM_Outputs"),
    variables=["prAdjust", "tasmaxAdjust", "tasminAdjust", "tasmeanAdjust"],
    frequencies=["Monthly", "Seasonal", "Annual"],

    scenarios=None,       # automatically discover available scenarios
    time_windows=None,    # automatically discover available time windows
    gcms=None,
    rcms=None,
    feature_ids=None,

    output_folder_name="Climate_Indicies",
    overwrite=False,
    plot_dpi=300,
    show_plots=False,
    verbose=True,
)

print(manifest["status"].value_counts(dropna=False))

The workflow preserves monthly/seasonal/annual time steps and produces
gridded NetCDF ensemble statistics, study-area CSV summaries, and
uncertainty plots.

6. Approximate climate-adjusted flood analysis

The flood workflow can run independently of the NARCliM THREDDS download
workflow.

Point spatial file

from pathlib import Path
from narclim_workflow import run_flood_hazard

OUTPUT_ROOT = Path(r"C:\NARCliM_Outputs")

TIME_WINDOWS = {
    "baseline":  (1985, 2014),
    "mid_term":  (2041, 2060),
    "long_term": (2081, 2100),
}

flood_results = run_flood_hazard(
    output_root=OUTPUT_ROOT,
    time_windows=TIME_WINDOWS,
    scenarios=["historical", "ssp245", "ssp370"],
    location_path=Path(r"C:\path\to\flood_locations.gpkg"),
    name_field="Name",
    context_field="Context",
    overwrite=False,
    verbose=True,
)

context_field should identify each point as Urban or Rural.
Default current reference events are 1, 10, and 100 year ARIs;
default durations are 1 hour for Urban and 24 hours for Rural. These
settings are configurable.

Principal outputs:

Flood_Hazard.gpkg
Flood_Hazard.csv
Flood_Runoff_Curves.csv

Downloaded/source information and QA files are retained under
NARCliM Data/Flood/.

Output structure

A typical project workspace is:

<OUTPUT_ROOT>/
├── NARCliM Data/
│   ├── prAdjust/
│   ├── tasmaxAdjust/
│   ├── tasminAdjust/
│   ├── TXge35/
│   ├── R99p/
│   ├── FFDIgt50/
│   ├── SPI12/
│   ├── sfcWindmax/              # optional retained storm source
│   ├── StormDaysGT89/           # threshold-dependent name
│   └── Flood/                   # BoM/ARR flood source + QA data
│
└── Climate_Indicies/
    ├── Temporal_Average/
    ├── Ensemble_Stats/
    ├── Hazards/
    │   ├── Ensemble_summary_all_variables.gpkg
    │   ├── Ensemble_summary_<variable>.gpkg
    │   ├── <dataset>_Hazards.gpkg
    │   ├── <dataset>_Hazards.csv
    │   └── Flood/
    │       ├── Flood_Hazard.gpkg
    │       ├── Flood_Hazard.csv
    │       └── Flood_Runoff_Curves.csv
    ├── Rain_Temp_Timeseries/
    │   └── <scenario>/<time_window>/<frequency>/
    │       ├── NetCDF/
    │       ├── CSV/
    │       └── Plots/
    └── Plot and Maps/

Climate_Indicies is retained because it is the current package
output-folder name. If the spelling is standardised to
Climate_Indices later, package configuration, examples, and
downstream paths should be updated together.

QA and interpretation

Temporal mean vs ensemble statistics

Raw climate time series
        ↓
Temporal mean within each GCM × RCM
        ↓
One value per grid cell per ensemble member
        ↓
Ensemble Min / [optional Mean] / Max across members

This processing order is important when interpreting outputs.

Polygon summaries

For polygon datasets, min, optional mean, and max are spatial
means of the corresponding ensemble rasters over each polygon. They are
not the minimum/maximum raster cells within that polygon.

Point summaries

For point datasets, ensemble statistics are sampled directly from the
raster cell containing each point.

Rainfall/temperature statistics

For gridded time-series outputs, ensemble statistics are calculated
pixel-by-pixel at each monthly, seasonal, or annual time step. For CSVs
and plots, each GCM × RCM member is first spatially averaged over the
study-area subset; ensemble statistics are then calculated across those
member-average series.

Flood interpretation

flood.py estimates climate-driven changes in equivalent event
frequency using rainfall-runoff matching. It should be interpreted as a
screening/comparative method and does not replace project-specific
hydrologic/hydraulic modelling.

Raster QA

Always inspect CRS, extent, resolution, valid-cell count, nodata
handling, and grid alignment as part of project QA.

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
confirming the licensing terms you want to use.

Citation

If you use or adapt this package in a project, please cite the package
and developer. A CITATION.cff file is included for GitHub citation
metadata.

Jabbar Khaledi
Data and Geospatial Analyst | Python Developer
Email: jabbarkhaledi88@gmail.com
Primary language: Python
Applications: NARCliM2.0 climate-data processing, climate-model
uncertainty analysis, rainfall/temperature time-series analysis,
storm-hazard analysis, approximate climate-adjusted flood-hazard
assessment, climate-risk assessment, and geospatial analysis.