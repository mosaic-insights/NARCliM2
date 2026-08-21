"""
NARCliM2 Workflow — Complete Example
====================================

Author / Developer
------------------
Jabbar Khaledi
Data and Geospatial Analyst | Python Developer

This file mirrors the content and workflow demonstrated in
NARCliM_Example.ipynb, but is formatted as a Python script suitable for a
GitHub repository.

PACKAGE OVERVIEW
----------------
The NARCliM2 Workflow package provides a flexible and repeatable workflow for
downloading, processing, summarising, and mapping NARCliM2 climate data from
NCI THREDDS.

The workflow operates through three main dimensions:

1. Temporal averaging
   For each variable × scenario × GCM × RCM × grid cell, values are averaged
   across all time steps in the selected horizon.

2. Ensemble uncertainty
   Temporal-average rasters from available GCM × RCM members are combined.
   Pixel-wise ensemble minimum, mean, and maximum are calculated.

3. Spatial summarisation
   Ensemble rasters are summarised for arbitrary polygon or point datasets.

For polygons, min/mean/max are the spatial means of the corresponding
ensemble-min/mean/max rasters. For points, values are sampled from the raster
cell containing the point.

StormDays
---------
StormDays is derived from daily NARCliM sfcWindmax.

The user controls the threshold with:

    storm_threshold_kmh=89.0

A threshold of 89 km/h creates StormDaysGT89. The original daily sfcWindmax
subset can also be retained with:

    save_storm_source=True

StormDaysGT* outputs use three decimal places because rare-event temporal and
ensemble averages can be smaller than 0.1 days/year. Most other climate values
are exported to one decimal place.

IMPORTANT
---------
Replace the example paths below with your own paths before running this file.
Large multi-variable daily-data runs can download/process substantial data.
"""

from pathlib import Path

from narclim_workflow import run_workflow
from narclim_workflow.uncertainty import run_uncertainty_workflow
from narclim_workflow.spatial_summary import run_spatial_summary
from narclim_workflow.ensemble_maps import plot_ensemble_maps


# =============================================================================
# USER PATHS
# =============================================================================

VECTOR_PATH = Path(r"C:\path\to\study_area.shp")
OUTPUT_ROOT = Path(r"C:\NARCliM_Outputs")

# Optional polygon/point dataset for spatial summaries.
FEATURE_PATH = Path(r"C:\path\to\planning_zones.shp")

CLIMATE_INDICES_ROOT = (
    OUTPUT_ROOT
    / "Climate_Indicies"
)


# =============================================================================
# USER-DEFINED TIME WINDOWS
#
# Historical NARCliM data generally cover 1951-2014.
# Projection scenarios generally cover 2015-2100.
# =============================================================================

TIME_WINDOWS = {
    "baseline": (1985, 2014),
    # "short_term": (2021, 2040),
    "mid_term": (2041, 2060),
    "long_term": (2081, 2100),
}


# =============================================================================
# WORKFLOW 1: DOWNLOAD / EXTRACT NARCliM DATA
#
# Select study area
#     ↓
# Select domain automatically
#     ↓
# Loop through variables, scenarios, models, and time windows
#     ↓
# Read NCI THREDDS datasets
#     ↓
# Subset time period
#     ↓
# Extract intersecting climate grid cells
#     ↓
# Save local NetCDF subsets
# =============================================================================

download_manifest = run_workflow(
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

    # "auto", "seaus_4km", or "aus_18"
    domain_key="auto",

    # None uses generated IDs such as ID1.
    id_field=None,

    # False is safer when re-running a large existing download.
    overwrite=False,

    print_traceback=True,

    # -----------------------------------------------------------------
    # Storm-specific settings.
    #
    # These are ignored for non-storm variables.
    # -----------------------------------------------------------------

    storm_threshold_kmh=89.0,

    # Retain the original daily sfcWindmax study-area subsets.
    save_storm_source=True,
)

print("\nDownload status summary:")
print(
    download_manifest[
        "status"
    ].value_counts(
        dropna=False
    )
)


# =============================================================================
# WORKFLOW 2: CLIMATE-MODEL UNCERTAINTY ANALYSIS
#
# Each GCM × RCM output
#         ↓
# Average each grid cell across all years/time steps in the selected horizon
#         ↓
# Save one temporal-average raster per ensemble member
#         ↓
# Combine all available ensemble-member rasters
#         ↓
# Calculate pixel-wise ensemble mean, minimum, and maximum
#         ↓
# Save three ensemble rasters
#         ↓
# Create pixel-centroid ensemble-summary GeoPackages
#
# IMPORTANT:
# If StormDays was downloaded using storm_threshold_kmh=89.0, the derived
# variable name used here is StormDaysGT89.
# =============================================================================

uncertainty_manifest = run_uncertainty_workflow(
    input_root=OUTPUT_ROOT,
    boundary_path=VECTOR_PATH,

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

print("\nUncertainty status summary:")
print(
    uncertainty_manifest[
        "status"
    ].value_counts(
        dropna=False
    )
)


# =============================================================================
# WORKFLOW 3: SPATIAL SUMMARY
#
# Read point or polygon dataset
#     ↓
# Discover available ensemble rasters
#     ↓
# Read climate/time metadata automatically
#     ↓
# Polygon?
#     ├── Yes → zonal spatial mean using touched cells for
#     │          ensemble min / mean / max
#     │
#     └── No → sample raster cell containing each point
#     ↓
# Attach min / mean / max to input features
#     ↓
# Add representative time_horizon_year
#     ↓
# Save one layer per climate variable in a GeoPackage
#     ↓
# Save a combined non-spatial CSV
#
# This can be used for planning zones, land use, districts, LGAs,
# catchments, subcatchments, sites, assets, or other user datasets.
# =============================================================================

spatial_outputs = run_spatial_summary(
    feature_path=FEATURE_PATH,
    climate_indices_root=CLIMATE_INDICES_ROOT,

    # Set to None to derive the output name from FEATURE_PATH.
    dataset_name="Planning_Zones",

    feature_id_field=None,

    # For polygons, include any raster cell touched by the polygon.
    all_touched=True,

    overwrite=True,
    verbose=True,
)

print("\nSpatial-summary outputs:")
for name, path in spatial_outputs.items():
    print(f"  {name}: {path}")


# ---------------------------------------------------------------------
# Optional restricted spatial-summary example.
# ---------------------------------------------------------------------
#
# restricted_outputs = run_spatial_summary(
#     feature_path=FEATURE_PATH,
#     climate_indices_root=CLIMATE_INDICES_ROOT,
#     variables=[
#         "TXge35",
#         "FFDIgt50",
#     ],
#     scenarios=[
#         "ssp245",
#         "ssp370",
#     ],
#     time_horizons=[
#         "mid_term",
#         "long_term",
#     ],
#     all_touched=True,
#     overwrite=True,
# )


# =============================================================================
# WORKFLOW 4: PLOTS AND MAPS
#
# Reads Ensemble_Stats GeoTIFFs, plots them using viridis, overlays the
# supplied study-area boundary, and saves PNG maps.
# =============================================================================

maps = plot_ensemble_maps(
    input_root=OUTPUT_ROOT,
    boundary_path=VECTOR_PATH,

    # None = all available variables/scenarios.
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

print(
    f"\nCreated/found {len(maps)} map(s)."
)
