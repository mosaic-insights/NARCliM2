# NARCliM2 Workflow

A reusable Python workflow for discovering, downloading, spatially subsetting, processing, summarising, and mapping **NARCliM2.0** climate data from the NCI THREDDS service.

The package is designed for climate-risk and geospatial applications where users need a repeatable workflow across study areas, climate variables, model ensembles, scenarios, and planning horizons.

## Author and Developer

**Jabbar Khaledi**  
*Data and Geospatial Analyst | Python Developer*

This package was designed and developed by **Jabbar Khaledi** to support reproducible spatial analysis of NARCliM2.0 climate projection data.

The workflow integrates NCI THREDDS data access, temporal processing, climate-model ensemble analysis, spatial summarisation, and mapping within a reusable Python package.

### Development scope

The package includes functionality for:

- automated access to NARCliM2.0 datasets through NCI THREDDS;
- automatic selection of an appropriate NARCliM spatial domain;
- spatial extraction for user-defined polygon or point study areas;
- processing of historical and future climate scenarios and time horizons;
- temporal averaging for individual GCM × RCM ensemble members;
- pixel-wise ensemble minimum, mean, and maximum calculations;
- derivation of configurable storm-wind threshold indicators;
- retention of original daily `sfcWindmax` storm-source data for QA and alternative analyses;
- raster, GeoPackage, CSV, and map outputs;
- spatial summarisation for planning zones, land use, administrative boundaries, catchments, sites, assets, and other spatial datasets.

### Design

**Reproducibility** — climate-data processing is implemented as repeatable Python workflows rather than manual GIS steps.

**Flexibility** — users can specify study areas, climate variables, scenarios, time horizons, GCMs, RCMs, and selected hazard parameters.

**Transparency** — intermediate products, temporal averages, ensemble statistics, metadata, and processing manifests are retained where useful.

**Spatial applicability** — climate-model outputs are converted to analysis-ready spatial products for planning, risk assessment, and GIS analysis.

---

## Package overview

The workflow processes NARCliM climate data through three main dimensions: **time**, **climate-model ensemble**, and **space**.

### 1. Temporal averaging — averaging through time

For each climate variable × scenario × GCM × RCM × grid cell, climate values are averaged across all time steps within the selected time horizon.

For example, a mid-term layer represented by 2041–2060 contains, for every ensemble member and grid cell, the average climate value over those years.

The result is **one temporal-average raster for each GCM × RCM ensemble member**.

### 2. Ensemble uncertainty — summarising across climate models

For each variable, scenario, and time horizon, the temporal-average rasters from all available GCM × RCM combinations are combined.

With 5 GCMs × 2 RCMs, each pixel can contain up to 10 ensemble-member values.

For every pixel independently, the package calculates:

- **Minimum (`min`)** — lowest value across ensemble members;
- **Mean (`mean`)** — average across ensemble members;
- **Maximum (`max`)** — highest value across ensemble members.

These produce three ensemble-statistic rasters representing the model range and central estimate at each grid cell.

Example filenames:

```text
TXge35_Min_historical_1985_2014.tif
TXge35_Mean_ssp245_2041_2060.tif
TXge35_Max_ssp370_2081_2100.tif
```

### 3. Spatial summarisation — applying climate information to features

The ensemble statistics can then be summarised for user-supplied spatial datasets such as:

- planning zones;
- land-use areas;
- districts;
- LGAs;
- catchments and subcatchments;
- properties;
- sites and assets;
- other polygon or point datasets.

For **polygon features**, zonal statistics are applied to the three ensemble rasters and the spatial mean of intersecting raster cells is calculated separately for ensemble minimum, mean, and maximum.

For **point features**, the minimum, mean, and maximum values are sampled from the raster cell containing each point.

The spatial-summary outputs also contain a representative `time_horizon_year`, for example:

| Period | Representative year |
|:--|--:|
| 1985–2014 | 2000 |
| 2041–2060 | 2050 |
| 2081–2100 | 2090 |

---

## Key features

- Download climate variables and indices directly from NCI THREDDS using OPeNDAP.
- Read `catalog.xml` rather than guessing remote filenames.
- Work with both **polygon** and **point** datasets.
- Process multiple input features automatically.
- Run multiple GCMs, RCMs, scenarios, variables, and time periods in a single workflow.
- Automatically clip/extract only climate cells relevant to the supplied study area.
- Automatically handle unavailable datasets by skipping them and reporting the reason.
- Calculate temporal-average rasters for each ensemble member.
- Calculate pixel-wise ensemble minimum, mean, and maximum rasters.
- Derive a configurable annual storm-days index from daily `sfcWindmax`.
- Save original daily storm-source subsets for QA if requested.
- Summarise ensemble hazard rasters for arbitrary polygon or point datasets.
- Save GeoPackage and CSV spatial-summary outputs.
- Produce `viridis` maps from ensemble rasters with the study-area boundary overlaid.
- Write processing/download manifests for QA and reproducibility.

---

## Installation

Clone the repository and install it in editable mode:

```bash
git clone <YOUR-REPOSITORY-URL>
cd narclim-workflow
python -m pip install -e .
```

For development with the example notebook:

```bash
python -m pip install -e ".[dev]"
```

Python **3.11+** is recommended.

---

## Domain selection

NARCliM2 provides data for two spatial domains used by this package.

| Domain | Approx. resolution | Coverage |
|:-------|:------------------:|:---------|
| **NARCliM2-0-SEAus-04** | ~4 km | South-East Australia |
| **AUS-18** | ~18 km | Australia |

The package supports three domain-selection modes:

| Option | Description |
|:------|:------------|
| `domain_key="auto"` | Automatically selects the appropriate domain based on the supplied study area. |
| `domain_key="seaus_4km"` | Forces the South-East Australia ~4 km domain. |
| `domain_key="aus_18"` | Forces the Australia-wide ~18 km domain. |

With `domain_key="auto"`:

- if all input features fall completely within the stored SEAus footprint, the workflow selects `NARCliM2-0-SEAus-04`;
- otherwise, it switches to `AUS-18`;
- the selected domain and reason are printed during execution.

---

## Variables configured in the package

| Variable | Description | NCI Output Branch | Processing Version | Frequency | AUS-18 | NARCliM2-0-SEAus-04 |
|:---------|:------------|:------------------|:-------------------|:---------:|:------:|:-------------------:|
| **prAdjust** | Daily bias-adjusted precipitation | `bias-adjusted-output` | `v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009` | `day` | ✓ Supported* | ✓ Supported |
| **tasmaxAdjust** | Daily bias-adjusted maximum temperature | `bias-adjusted-output` | `v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009` | `day` | ✓ Supported* | ✓ Supported |
| **tasminAdjust** | Daily bias-adjusted minimum temperature | `bias-adjusted-output` | `v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009` | `day` | ✓ Supported* | ✓ Supported |
| **TXge35** | Yearly number of days with maximum temperature ≥ 35°C | `bias-adjusted-output` | `v1-r1` | `yr` | ✓ Supported* | ✓ Supported |
| **TNlt2** | Yearly number of days with minimum temperature < 2°C | `bias-adjusted-output` | `v1-r1` | `yr` | ✓ Supported* | ✓ Supported |
| **FFDIgt50** | Yearly number of days with FFDI ≥ 50 | `DD` | `v1-r1` | `yr` | ✓ Supported* | ✓ Supported |
| **R20mm** | Yearly number of days with precipitation ≥ 20 mm | `DD` | `v1-r1` | `yr` | ✓ Supported* | ✓ Supported |
| **R99p** | Yearly total precipitation from extremely wet days | `DD` | `v1-r1` | `yr` | ✓ Supported* | ✓ Supported |
| **SPI12** | 12-month Standardised Precipitation Index | `DD` | `v1-r1` | `mon` | ✓ Supported* | ✓ Supported |
| **StormDays** | Derived annual count of days where daily maximum near-surface wind exceeds a user-defined threshold | `DD` (`sfcWindmax`) | `v1-r1` | source: `day`; derived: `yr` | ✓ Supported* | ✓ Supported |

> **Supported*** means the package can construct the expected NCI THREDDS catalogue path for the selected domain. Actual availability still depends on the requested domain, scenario, GCM, RCM, and period. Unavailable datasets are skipped and reported.

### Storm processing

`StormDays` is a **derived climate-hazard indicator**, not a directly downloaded annual NARCliM index.

The workflow reads daily NARCliM **maximum near-surface wind speed**, `sfcWindmax`, in `m s-1`.

The source data are spatially subset to the supplied study area. If:

```python
save_storm_source=True
```

the original daily subset is retained under:

```text
NARCliM Data/
└── sfcWindmax/
    ├── historical/
    ├── ssp245/
    └── ssp370/
```

The user specifies the threshold in km/h:

```python
storm_threshold_kmh=89.0
```

The package converts the threshold to m/s and calculates, for every grid cell and model year:

> **number of days where daily `sfcWindmax` > threshold**

A threshold of 89 km/h creates:

```text
StormDaysGT89
```

and stores the annual derived data under:

```text
NARCliM Data/
└── StormDaysGT89/
```

A different threshold changes the derived output name automatically. For example:

```python
storm_threshold_kmh=63.0
```

produces:

```text
StormDaysGT63
```

Retaining `sfcWindmax` is useful for QA and for testing alternative thresholds without needing to download the same daily wind source again.

### Storm precision

High-wind exceedances can be rare. After temporal and ensemble averaging, meaningful storm frequencies can be small, for example:

```text
0.033 days/year
```

To avoid converting these values to `0.0`, `StormDaysGT*` products retain **three decimal places**, while most other exported climate values use **one decimal place**.

---

## Flexible variable selection

Any combination of configured variables can be processed in one run:

```python
variables=[
    "prAdjust",
    "tasmaxAdjust",
    "tasminAdjust",
    "TXge35",
    "R99p",
    "FFDIgt50",
    "SPI12",
    "StormDays",
]
```

Only requested variables are processed.

---

## Climate-model selection

Users can choose any combination of GCMs:

```python
gcms=[
    "ACCESS-ESM1-5",
    "EC-Earth3-Veg",
    "MPI-ESM1-2-HR",
    "NorESM2-MM",
    "UKESM1-0-LL",
]
```

and either or both RCMs:

```python
rcms=[
    "NARCliM2-0-WRF412R3",
    "NARCliM2-0-WRF412R5",
]
```

The workflow loops through the requested combinations automatically.

---

## Scenario selection

Any configured combination can be requested:

```python
scenarios=[
    "historical",
    "ssp126",
    "ssp245",
    "ssp370",
]
```

Variables unavailable for a scenario are skipped without stopping the full batch. For example, the package configuration treats `SPI12` as unavailable for `historical`.

---

## Flexible time windows

Users define planning horizons with a dictionary:

```python
TIME_WINDOWS = {
    "baseline": (1985, 2014),
    "short_term": (2021, 2040),
    "mid_term": (2041, 2060),
    "long_term": (2081, 2100),
}
```

Historical data are generally available within **1951–2014**, while projection scenarios cover **2015–2100**. The requested windows should be confirmed for the intended analysis.

---

## Spatial flexibility

The download workflow accepts:

- polygon layers such as catchments, LGAs, planning zones, properties, and study boundaries;
- point layers such as monitoring sites, assets, and infrastructure.

If multiple features are supplied, each feature is processed independently and stored using a feature ID.

---

## Automatic handling of unavailable data

Rather than terminating a large batch run, the workflow:

- checks the requested dataset/catalogue;
- skips unavailable datasets;
- reports the reason;
- continues with remaining combinations.

This is particularly useful for runs involving many variables, scenarios, GCMs, and RCMs.

---

## Overall workflow

```text
workflow.py
    ↓
Download/extract NARCliM data
    ↓
uncertainty.py
    ↓
Temporal mean for each model
    ↓
Pixel-wise ensemble Min / Mean / Max
    ↓
spatial_summary.py
    ↓
Summarise hazard rasters for arbitrary polygon or point datasets
    ↓
ensemble_maps.py
    ↓
Create maps from Ensemble_Stats rasters
```

### Download workflow

```text
Select study area
    ↓
Select domain automatically or explicitly
    ↓
Loop through variables, scenarios, models, and time windows
    ↓
Read NCI THREDDS catalogues
    ↓
Open remote NetCDF via OPeNDAP
    ↓
Subset requested time period
    ↓
Extract intersecting climate grid cells
    ↓
Save local NetCDF subsets
```

### Uncertainty workflow

```text
Each GCM × RCM output
        ↓
Average each grid cell across all years/time steps in the selected horizon
        ↓
Save one temporal-average raster per ensemble member
        ↓
Combine all available ensemble-member rasters
        ↓
Calculate pixel-wise ensemble mean, minimum, and maximum
        ↓
Save three ensemble rasters
        ↓
Create pixel-centroid ensemble-summary GeoPackages
```

### Spatial-summary workflow

```text
Read polygon or point dataset
    ↓
Discover available ensemble rasters
    ↓
Read climate/time metadata automatically
    ↓
Polygon?
    ├── Yes → zonal spatial mean using touched raster cells
    │          for ensemble min / mean / max
    │
    └── No → sample raster cell containing each point
    ↓
Attach min / mean / max to input features
    ↓
Add representative time_horizon_year
    ↓
Save one layer per climate variable in a GeoPackage
    ↓
Save a combined non-spatial CSV
```

---

## Output structure

A typical workspace contains:

```text
<OUTPUT_ROOT>/
│
├── NARCliM Data/
│   ├── prAdjust/
│   ├── tasmaxAdjust/
│   ├── tasminAdjust/
│   ├── TXge35/
│   ├── R99p/
│   ├── FFDIgt50/
│   ├── SPI12/
│   ├── sfcWindmax/          # optional retained daily storm source
│   └── StormDaysGT89/       # name changes with threshold
│
└── Climate_Indicies/
    ├── Temporal_Average/
    ├── Ensemble_Stats/
    ├── Hazards/
    │   ├── Ensemble_summary_all_variables.gpkg
    │   ├── Ensemble_summary_<variable>.gpkg
    │   ├── <dataset>_Hazards.gpkg
    │   └── <dataset>_Hazards.csv
    └── Plot and Maps/
```

> `Climate_Indicies` is retained here because it is the current package output-folder name. If the package later standardises the spelling to `Climate_Indices`, update the examples and downstream paths together.

---

# Examples

The following examples mirror the complete workflow shown in `NARCliM_Example.ipynb`. Replace the example paths with your own project paths.

## 1. Download NARCliM data

```python
from pathlib import Path

from narclim_workflow import run_workflow


VECTOR_PATH = Path(r"C:\path\to\study_area.shp")
OUTPUT_ROOT = Path(r"C:\NARCliM_Outputs")

TIME_WINDOWS = {
    "baseline": (1985, 2014),
    # "short_term": (2021, 2040),
    "mid_term": (2041, 2060),
    "long_term": (2081, 2100),
}

manifest = run_workflow(
    vector_path=VECTOR_PATH,
    output_root=OUTPUT_ROOT,
    time_windows=TIME_WINDOWS,
    variables=[
        "prAdjust",
        "tasmaxAdjust",
        "tasminAdjust",
        "TXge35",
        "R99p",
        "FFDIgt50",
        "SPI12",
        "StormDays",
    ],
    gcms=[
        "ACCESS-ESM1-5",
        "EC-Earth3-Veg",
        "MPI-ESM1-2-HR",
        "NorESM2-MM",
        "UKESM1-0-LL",
    ],
    scenarios=[
        "historical",
        "ssp245",
        "ssp370",
    ],
    rcms=[
        "NARCliM2-0-WRF412R3",
        "NARCliM2-0-WRF412R5",
    ],
    domain_key="auto",
    id_field=None,
    overwrite=False,
    print_traceback=True,

    # Storm-specific settings. Ignored for non-storm variables.
    storm_threshold_kmh=89.0,
    save_storm_source=True,
)

print(manifest["status"].value_counts(dropna=False))
```

> Large daily-variable runs can be substantial. If existing outputs should be reused, keep `overwrite=False`.

## 2. Climate-model uncertainty analysis

For `StormDays`, use the dynamically generated variable name corresponding to the threshold used during download. With `storm_threshold_kmh=89.0`, that variable is `StormDaysGT89`.

```python
from pathlib import Path

from narclim_workflow.uncertainty import run_uncertainty_workflow


INPUT_ROOT = Path(r"C:\NARCliM_Outputs")
BOUNDARY_PATH = Path(r"C:\path\to\study_area.shp")

manifest = run_uncertainty_workflow(
    input_root=INPUT_ROOT,
    boundary_path=BOUNDARY_PATH,
    variables=[
        "TXge35",
        "R20mm",
        "FFDIgt50",
        "SPI12",
        "StormDaysGT89",
    ],
    scenarios=[
        "historical",
        "ssp126",
        "ssp245",
        "ssp370",
    ],
    time_windows=[
        "baseline",
        "short_term",
        "mid_term",
        "long_term",
    ],
    gcms=[
        "ACCESS-ESM1-5",
        "EC-Earth3-Veg",
        "MPI-ESM1-2-HR",
        "NorESM2-MM",
        "UKESM1-0-LL",
    ],
    rcms=[
        "NARCliM2-0-WRF412R3",
        "NARCliM2-0-WRF412R5",
    ],
    feature_ids=None,
    boundary_id_field=None,
    output_folder_name="Climate_Indicies",
    overwrite=False,
    verbose=True,
)

print("\nStatus summary:")
print(manifest["status"].value_counts(dropna=False))
```

The uncertainty workflow first calculates the temporal average within each GCM × RCM member and then calculates ensemble min/mean/max **across model-member temporal averages**. Therefore ensemble min/max are model-ensemble bounds, not the minimum/maximum individual year within the horizon.

## 3. Spatial summary for planning zones or other datasets

```python
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
    overwrite=True,
    verbose=True,
)
```

To let the module derive the output name from the input dataset, use:

```python
dataset_name=None
```

You can also restrict processing:

```python
outputs = run_spatial_summary(
    feature_path=FEATURE_PATH,
    climate_indices_root=CLIMATE_INDICES_ROOT,
    variables=["TXge35", "FFDIgt50"],
    scenarios=["ssp245", "ssp370"],
    time_horizons=["mid_term", "long_term"],
    all_touched=True,
    overwrite=True,
)
```

For a `Planning_Zones` dataset, outputs include:

```text
Planning_Zones_Hazards.gpkg
Planning_Zones_Hazards.csv
```

Climate values use one decimal place for most variables and three decimals for `StormDaysGT*`.

## 4. Plot ensemble rasters

```python
from pathlib import Path

from narclim_workflow.ensemble_maps import plot_ensemble_maps


INPUT_ROOT = Path(r"C:\NARCliM_Outputs")
BOUNDARY_PATH = Path(r"C:\path\to\study_area.shp")

maps = plot_ensemble_maps(
    input_root=INPUT_ROOT,
    boundary_path=BOUNDARY_PATH,
    variables=None,
    scenarios=None,
    statistics=[
        "min",
        "mean",
        "max",
    ],
    cmap="viridis",
    dpi=300,
    overwrite=True,
    show=False,
    verbose=True,
)

print(f"Created/found {len(maps)} map(s).")
```

Maps are saved under:

```text
Climate_Indicies/
└── Plot and Maps/
```

The study-area boundary is reprojected to the raster CRS before plotting.

---

## QA and interpretation notes

### Temporal mean vs ensemble statistics

The processing order is:

```text
raw time series
    ↓
temporal mean within each GCM × RCM
    ↓
one value per pixel per ensemble member
    ↓
ensemble minimum / mean / maximum across members
```

This distinction is important when interpreting the output.

### Polygon summaries

For polygon datasets, the output `min`, `mean`, and `max` fields are spatial means of the corresponding ensemble rasters over the polygon. They are not the minimum/maximum raster cell within the polygon.

### Point summaries

For points, `min`, `mean`, and `max` are sampled directly from the intersecting ensemble raster cell.

### Raster resolution and CRS

The package retains the NARCliM grid information and writes ensemble raster products in a GIS-compatible form. Always inspect CRS, extent, resolution, valid-cell count, and alignment as part of project QA.

---

## Repository structure

A recommended GitHub layout is:

```text
narclim-workflow/
├── README.md
├── pyproject.toml
├── example.py
├── NARCliM_Example.ipynb
├── CITATION.cff
├── CHANGELOG.md
├── .gitignore
├── src/
│   └── narclim_workflow/
│       ├── __init__.py
│       ├── catalog.py
│       ├── config.py
│       ├── domain.py
│       ├── spatial.py
│       ├── workflow.py
│       ├── uncertainty.py
│       ├── spatial_summary.py
│       └── ensemble_maps.py
└── tests/
```

For a public release, also choose and add an appropriate open-source `LICENSE` after confirming the licensing terms you want to use.

---

## Citation

If you use or adapt this package in a project, please cite the package and developer. A `CITATION.cff` file is included to make GitHub citation metadata available.

**Developer:** Jabbar Khaledi  (Email: jabbarkhaledi88@gmail.com) 
**Primary language:** Python  
**Application:** NARCliM2.0 climate-data processing, climate-model uncertainty analysis, climate-risk assessment, and geospatial analysis
