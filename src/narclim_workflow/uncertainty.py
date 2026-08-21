from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import traceback
from typing import Iterable
from .config import VARIABLES
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
import xarray as xr
from pyproj import CRS as PyprojCRS
from rasterio.crs import CRS as RasterioCRS
from rasterio.features import rasterize
from shapely.geometry import Point, Polygon

UNCERTAINTY_CODE_VERSION = "2026-08-21-storm-precision-v7"

@dataclass(frozen=True)
class ModelLayer:
    variable: str
    scenario: str
    window: str
    gcm: str
    rcm: str
    feature_id: str
    source_path: Path
    data: xr.DataArray


_STORM_DERIVED_RE = re.compile(
    r"^StormDaysGT\d+(?:p\d+)?$",
    flags=re.IGNORECASE,
)


def _is_storm_derived_variable(
    variable: str | None,
) -> bool:
    """Return True for dynamic storm-day variables such as StormDaysGT89."""

    if variable is None:
        return False

    return bool(
        _STORM_DERIVED_RE.match(
            str(variable)
        )
    )


def _output_decimals(
    variable: str | None,
) -> int:
    """
    Return export precision for a climate variable.

    Storm-day frequencies use 3 decimals because temporal/ensemble averages
    can be very small but meaningful (for example 1/30 = 0.033 days/year).
    Other climate variables keep the one-decimal package convention.
    """

    if _is_storm_derived_variable(
        variable
    ):
        return 3

    return 1



def _centres_to_corners(
    values: np.ndarray,
) -> np.ndarray:
    """
    Convert a 2-D array of climate-cell centre coordinates into
    approximate cell-corner coordinates.
    """

    values = np.asarray(
        values,
        dtype=np.float64,
    )

    if values.ndim != 2:
        raise ValueError(
            "Cell-centre coordinates must be two-dimensional."
        )

    n_rows, n_cols = values.shape

    if n_rows < 2 or n_cols < 2:
        raise ValueError(
            "At least a 2 x 2 grid is required."
        )

    corners = np.empty(
        (
            n_rows + 1,
            n_cols + 1,
        ),
        dtype=np.float64,
    )

    # Interior corners.
    corners[
        1:n_rows,
        1:n_cols,
    ] = (
        values[:-1, :-1]
        + values[:-1, 1:]
        + values[1:, :-1]
        + values[1:, 1:]
    ) / 4.0

    # Top edge.
    top_midpoints = (
        values[0, :-1]
        + values[0, 1:]
    ) / 2.0

    corners[
        0,
        1:n_cols,
    ] = (
        2.0 * top_midpoints
        - corners[1, 1:n_cols]
    )

    # Bottom edge.
    bottom_midpoints = (
        values[-1, :-1]
        + values[-1, 1:]
    ) / 2.0

    corners[
        n_rows,
        1:n_cols,
    ] = (
        2.0 * bottom_midpoints
        - corners[
            n_rows - 1,
            1:n_cols,
        ]
    )

    # Left edge.
    left_midpoints = (
        values[:-1, 0]
        + values[1:, 0]
    ) / 2.0

    corners[
        1:n_rows,
        0,
    ] = (
        2.0 * left_midpoints
        - corners[
            1:n_rows,
            1,
        ]
    )

    # Right edge.
    right_midpoints = (
        values[:-1, -1]
        + values[1:, -1]
    ) / 2.0

    corners[
        1:n_rows,
        n_cols,
    ] = (
        2.0 * right_midpoints
        - corners[
            1:n_rows,
            n_cols - 1,
        ]
    )

    # Outer corners.
    corners[0, 0] = (
        corners[0, 1]
        + corners[1, 0]
        - corners[1, 1]
    )

    corners[0, n_cols] = (
        corners[0, n_cols - 1]
        + corners[1, n_cols]
        - corners[1, n_cols - 1]
    )

    corners[n_rows, 0] = (
        corners[n_rows - 1, 0]
        + corners[n_rows, 1]
        - corners[n_rows - 1, 1]
    )

    corners[n_rows, n_cols] = (
        corners[n_rows - 1, n_cols]
        + corners[n_rows, n_cols - 1]
        - corners[n_rows - 1, n_cols - 1]
    )

    return corners

def _native_output_resolution(
    source_ds: xr.Dataset,
) -> tuple[float, float]:
    """
    Return the package output pixel size in degrees.

    For SEAus-04 the package preserves the requested 0.04 x 0.04 degree
    grid in the GIS-ready EPSG:4326 outputs.

    For other NARCliM domains, estimate the native rotated-grid spacing
    from rlon/rlat and round to two decimal places, matching the coordinate
    precision used in the downloaded NARCliM subsets. This avoids deriving
    output resolution from geographic lon/lat distortion.
    """

    domain_id = str(
        source_ds.attrs.get(
            "domain_id",
            "",
        )
    )

    if domain_id == "NARCliM2-0-SEAus-04":
        return 0.04, 0.04

    if "rlon" not in source_ds.coords or "rlat" not in source_ds.coords:
        raise ValueError(
            f"Cannot determine native output resolution for domain '{domain_id}': "
            "rlon/rlat coordinates are missing."
        )

    rlon = np.asarray(source_ds["rlon"].values, dtype=np.float64)
    rlat = np.asarray(source_ds["rlat"].values, dtype=np.float64)

    if rlon.ndim != 1 or rlat.ndim != 1 or rlon.size < 2 or rlat.size < 2:
        raise ValueError(
            f"Cannot determine native output resolution for domain '{domain_id}'."
        )

    # Fit the nominal native spacing, then round to the same two-decimal
    # precision used by the saved NARCliM rotated coordinates.
    x_index = np.arange(rlon.size, dtype=np.float64)
    y_index = np.arange(rlat.size, dtype=np.float64)

    x_slope, _ = np.polyfit(x_index, rlon, 1)
    y_slope, _ = np.polyfit(y_index, rlat, 1)

    x_resolution = round(abs(float(x_slope)), 2)
    y_resolution = round(abs(float(y_slope)), 2)

    if x_resolution <= 0 or y_resolution <= 0:
        raise ValueError(
            f"Invalid native output resolution for domain '{domain_id}'."
        )

    return x_resolution, y_resolution


def _rasterize_curvilinear_to_wgs84(
    *,
    data: xr.DataArray,
    source_ds: xr.Dataset,
    nodata: float = -9999.0,
) -> tuple[np.ndarray, object]:
    """
    Convert selected NARCliM climate cells into a GIS-ready EPSG:4326
    raster while preserving the configured/native output resolution.

    Spatial logic
    -------------
    * NCI's 2-D lon/lat coordinates define the true geographic positions.
    * Each finite source climate cell is converted to a geographic polygon.
    * The complete source-cell polygon is rasterised with all_touched=True.
    * SEAus-04 outputs use exactly 0.04 x 0.04 degree pixels.
    * The output extent is snapped outward to that fixed grid so no selected
      climate-cell footprint is lost.

    Note that output width/height can differ from the source array shape.
    That is necessary because a rotated source grid cannot, in general, keep
    the same extent, same dimensions, and a fixed square EPSG:4326 pixel size
    simultaneously.
    """

    if "lon" not in source_ds.variables:
        raise ValueError(
            "Source dataset does not contain 2-D lon coordinates."
        )

    if "lat" not in source_ds.variables:
        raise ValueError(
            "Source dataset does not contain 2-D lat coordinates."
        )

    data_2d = data.squeeze(drop=True)

    y_dim, x_dim = _spatial_dims(data_2d)
    data_2d = data_2d.transpose(y_dim, x_dim)

    values = np.asarray(
        data_2d.values,
        dtype=np.float64,
    )

    lon = np.asarray(
        source_ds["lon"].values,
        dtype=np.float64,
    )

    lat = np.asarray(
        source_ds["lat"].values,
        dtype=np.float64,
    )

    if lon.shape != values.shape:
        raise ValueError(
            f"lon shape {lon.shape} does not match data shape {values.shape}."
        )

    if lat.shape != values.shape:
        raise ValueError(
            f"lat shape {lat.shape} does not match data shape {values.shape}."
        )

    n_rows, n_cols = values.shape

    # Build geographic cell corners from NCI's true 2-D cell centres.
    lon_corners = _centres_to_corners(lon)
    lat_corners = _centres_to_corners(lat)

    shapes = []
    selected_source_cells = 0

    west_values = []
    south_values = []
    east_values = []
    north_values = []

    for row in range(n_rows):
        for col in range(n_cols):
            value = values[row, col]

            # NaN means this native climate cell was outside the selected
            # study geometry and should not be transferred to the output.
            if not np.isfinite(value):
                continue

            coordinates = [
                (
                    lon_corners[row, col],
                    lat_corners[row, col],
                ),
                (
                    lon_corners[row, col + 1],
                    lat_corners[row, col + 1],
                ),
                (
                    lon_corners[row + 1, col + 1],
                    lat_corners[row + 1, col + 1],
                ),
                (
                    lon_corners[row + 1, col],
                    lat_corners[row + 1, col],
                ),
            ]

            if not all(
                np.isfinite(x) and np.isfinite(y)
                for x, y in coordinates
            ):
                continue

            cell_polygon = Polygon(coordinates)

            if not cell_polygon.is_valid:
                cell_polygon = cell_polygon.buffer(0)

            if cell_polygon.is_empty:
                continue

            shapes.append(
                (
                    cell_polygon,
                    float(value),
                )
            )

            selected_source_cells += 1

            minx, miny, maxx, maxy = cell_polygon.bounds
            west_values.append(minx)
            south_values.append(miny)
            east_values.append(maxx)
            north_values.append(maxy)

    if not shapes:
        raise ValueError(
            "No valid selected NARCliM climate cells were available "
            "for rasterisation."
        )

    # Fixed/native output pixel size. For SEAus this is exactly 0.04 x 0.04.
    x_resolution, y_resolution = _native_output_resolution(source_ds)

    # Use the complete selected-cell footprint, then snap the extent OUTWARD
    # to the fixed output grid so edge cells are not clipped away.
    raw_west = float(min(west_values))
    raw_south = float(min(south_values))
    raw_east = float(max(east_values))
    raw_north = float(max(north_values))

    west = float(
        np.floor(raw_west / x_resolution)
        * x_resolution
    )

    east = float(
        np.ceil(raw_east / x_resolution)
        * x_resolution
    )

    south = float(
        np.floor(raw_south / y_resolution)
        * y_resolution
    )

    north = float(
        np.ceil(raw_north / y_resolution)
        * y_resolution
    )

    width = max(
        1,
        int(
            round(
                (east - west)
                / x_resolution
            )
        ),
    )

    height = max(
        1,
        int(
            round(
                (north - south)
                / y_resolution
            )
        ),
    )

    transform = from_origin(
        west,
        north,
        x_resolution,
        y_resolution,
    )

    output = rasterize(
        shapes=shapes,
        out_shape=(
            height,
            width,
        ),
        transform=transform,
        fill=nodata,
        all_touched=True,
        dtype="float32",
    )

    valid_output_pixels = int(
        np.count_nonzero(
            output != nodata
        )
    )

    domain_id = str(
        source_ds.attrs.get(
            "domain_id",
            "unknown",
        )
    )

    return output, transform

def _log(message: str, *, verbose: bool) -> None:
    if verbose:
        print(message, flush=True)


def _safe_name(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r'[<>:"/\\|?*]+', "_", text)
    text = re.sub(r"\s+", "_", text)
    return text.strip("._") or "unnamed"


def _find_time_name(data: xr.DataArray) -> str | None:
    for name in list(data.dims) + list(data.coords):
        if name.lower() == "time":
            return name
    return None


def _choose_data_variable(ds: xr.Dataset, preferred: str) -> str:
    for name in ds.data_vars:
        if name.lower() == preferred.lower():
            return name

    candidates = [
        name
        for name in ds.data_vars
        if "bound" not in name.lower()
        and not name.lower().endswith("_bnds")
        and name.lower() != "crs"
    ]

    if len(candidates) == 1:
        return candidates[0]

    raise KeyError(
        f"Could not uniquely identify variable '{preferred}'. "
        f"Available data variables: {list(ds.data_vars)}"
    )


def _spatial_dims(data: xr.DataArray) -> tuple[str, str]:
    dim_lookup = {dim.lower(): dim for dim in data.dims}

    for y_name, x_name in [
        ("rlat", "rlon"),
        ("lat", "lon"),
        ("latitude", "longitude"),
        ("y", "x"),
    ]:
        if y_name in dim_lookup and x_name in dim_lookup:
            return dim_lookup[y_name], dim_lookup[x_name]

    non_time_dims = [
        dim
        for dim in data.dims
        if dim.lower() not in {"time", "bnds", "bounds", "crs"}
        and data.sizes.get(dim, 1) > 1
    ]

    if len(non_time_dims) == 2:
        return non_time_dims[0], non_time_dims[1]

    raise ValueError(f"Could not identify horizontal dimensions from {data.dims}")


def _crs_from_dataset(ds: xr.Dataset, data: xr.DataArray) -> CRS | None:
    """Try to construct a CRS from a dataset's CF grid-mapping metadata."""

    grid_mapping_name = data.attrs.get("grid_mapping")
    candidate_names = []

    if grid_mapping_name:
        candidate_names.append(grid_mapping_name)

    if "crs" not in candidate_names:
        candidate_names.append("crs")

    for variable_name in candidate_names:
        if variable_name not in ds.variables:
            continue

        attrs = dict(ds[variable_name].attrs)

        if attrs.get("spatial_ref"):
            return CRS.from_wkt(attrs["spatial_ref"])

        if attrs.get("crs_wkt"):
            return CRS.from_wkt(attrs["crs_wkt"])

        try:
            return CRS.from_cf(attrs)
        except Exception:
            continue

    return None

def _get_narclim_fallback_crs(
    source_ds: xr.Dataset,
) -> RasterioCRS | None:
    """
    Reconstruct the standard NARCliM2 rotated-pole CRS when an NCI
    NetCDF file does not contain a CF grid-mapping variable.

    This is required because some NCI files, for example some
    FFDIgt50 SSP245 datasets, contain rlon/rlat and lon/lat coordinates
    but do not contain a 'crs' variable or a grid_mapping attribute.

    The fallback is only used when the dataset can be confidently
    identified as a NARCliM rotated-pole dataset.
    """

    # ------------------------------------------------------------------
    # The rotated-grid coordinates must exist.
    # ------------------------------------------------------------------

    required_coordinates = {
        "rlon",
        "rlat",
        "lon",
        "lat",
    }

    if not required_coordinates.issubset(
        source_ds.variables
    ):
        return None

    # ------------------------------------------------------------------
    # Check that this really is a supported NARCliM domain.
    # ------------------------------------------------------------------

    domain_id = str(
        source_ds.attrs.get(
            "domain_id",
            "",
        )
    )

    grid_description = str(
        source_ds.attrs.get(
            "grid",
            "",
        )
    )

    supported_domains = {
        "NARCliM2-0-SEAus-04",
        "AUS-18",
    }

    if domain_id not in supported_domains:
        return None

    if (
        "rotated-pole" not in grid_description.lower()
        and
        "rotated pole" not in grid_description.lower()
    ):
        return None

    # ------------------------------------------------------------------
    # Standard NARCliM2 rotated-pole definition.
    #
    # These values were confirmed from the working NCI NetCDF files.
    # ------------------------------------------------------------------

    cf_mapping = {
        "grid_mapping_name":
            "rotated_latitude_longitude",

        "grid_north_pole_latitude":
            60.31,

        "grid_north_pole_longitude":
            141.38,

        "north_pole_grid_longitude":
            0.0,

        "earth_radius":
            6370000.0,
    }

    pyproj_crs = PyprojCRS.from_cf(
        cf_mapping
    )

    return RasterioCRS.from_wkt(
        pyproj_crs.to_wkt()
    )

def _get_crs(
    source_ds: xr.Dataset,
    data: xr.DataArray,
) -> RasterioCRS:
    """
    Obtain a usable CRS for a NARCliM raster.

    Search order
    ------------
    1. Climate-variable grid_mapping metadata.
    2. Local 'crs' variable.
    3. Original NCI OPeNDAP source.
    4. Known NARCliM rotated-pole CRS fallback.

    The fourth method is required for NCI datasets that genuinely lack
    CF CRS metadata.
    """

    # ==================================================================
    # 1. GRID-MAPPING REFERENCED BY THE CLIMATE VARIABLE
    # ==================================================================

    grid_mapping_name = data.attrs.get(
        "grid_mapping"
    )

    if (
        grid_mapping_name
        and grid_mapping_name in source_ds.variables
    ):
        try:
            cf_attrs = dict(
                source_ds[
                    grid_mapping_name
                ].attrs
            )

            pyproj_crs = PyprojCRS.from_cf(
                cf_attrs
            )

            return RasterioCRS.from_wkt(
                pyproj_crs.to_wkt()
            )

        except Exception:
            pass

    # ==================================================================
    # 2. LOCAL CRS VARIABLE
    # ==================================================================

    if "crs" in source_ds.variables:
        try:
            cf_attrs = dict(
                source_ds["crs"].attrs
            )

            pyproj_crs = PyprojCRS.from_cf(
                cf_attrs
            )

            return RasterioCRS.from_wkt(
                pyproj_crs.to_wkt()
            )

        except Exception:
            pass

    # ==================================================================
    # 3. ORIGINAL NCI DATASET
    # ==================================================================

    source_opendap = source_ds.attrs.get(
        "source_opendap"
    )

    if source_opendap:

        try:
            with xr.open_dataset(
                source_opendap,
                engine="netcdf4",
                decode_times=False,
                chunks=None,
                cache=False,
            ) as remote_ds:

                # Try the same variable name first.
                if data.name in remote_ds.data_vars:
                    remote_var = remote_ds[
                        data.name
                    ]

                    remote_mapping = (
                        remote_var.attrs.get(
                            "grid_mapping"
                        )
                    )

                    if (
                        remote_mapping
                        and remote_mapping
                        in remote_ds.variables
                    ):
                        cf_attrs = dict(
                            remote_ds[
                                remote_mapping
                            ].attrs
                        )

                        pyproj_crs = (
                            PyprojCRS.from_cf(
                                cf_attrs
                            )
                        )

                        return (
                            RasterioCRS.from_wkt(
                                pyproj_crs.to_wkt()
                            )
                        )

                # Also try a generic CRS variable.
                if "crs" in remote_ds.variables:

                    cf_attrs = dict(
                        remote_ds["crs"].attrs
                    )

                    pyproj_crs = (
                        PyprojCRS.from_cf(
                            cf_attrs
                        )
                    )

                    return RasterioCRS.from_wkt(
                        pyproj_crs.to_wkt()
                    )

        except Exception:
            # Remote access is a fallback only.
            # Do not terminate the workflow here.
            pass

    # ==================================================================
    # 4. NARCLIM FALLBACK
    # ==================================================================

    fallback_crs = (
        _get_narclim_fallback_crs(
            source_ds
        )
    )

    if fallback_crs is not None:

        print(
            "[CRS FALLBACK] Source NetCDF has no usable CF CRS metadata; "
            "using the standard NARCliM2 rotated-pole CRS.",
            flush=True,
        )

        return fallback_crs

    # ==================================================================
    # NOTHING WORKED
    # ==================================================================

    raise ValueError(
        "No usable CF grid-mapping metadata was found and the "
        "dataset could not be identified as a supported NARCliM "
        "rotated-pole grid."
    )

def _fit_regular_axis(
    values: np.ndarray,
    *,
    axis_name: str,
) -> tuple[np.ndarray, float, float]:
    """
    Reconstruct a regular coordinate axis from rounded NARCliM coordinates.

    NARCliM rotated-grid coordinates may be stored with only two decimal
    places. Although the underlying grid is regular, this rounding can make
    adjacent coordinate differences alternate between values such as 0.03
    and 0.04.

    A linear model is fitted to the coordinate values:

        coordinate = intercept + grid_index * spacing

    Parameters
    ----------
    values
        One-dimensional coordinate values.

    axis_name
        Coordinate name used in validation and error messages.

    Returns
    -------
    fitted_values
        Reconstructed regular cell-centre coordinates.

    resolution
        Absolute fitted grid spacing.

    maximum_residual
        Maximum absolute difference between the stored and fitted values.
    """

    values = np.asarray(
        values,
        dtype=np.float64,
    )

    if values.ndim != 1:
        raise ValueError(
            f"{axis_name} must be one-dimensional. "
            f"Received shape {values.shape}."
        )

    if values.size < 2:
        raise ValueError(
            f"{axis_name} requires at least two coordinate values."
        )

    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"{axis_name} contains non-finite coordinate values."
        )

    differences = np.diff(values)

    increasing = bool(
        np.all(differences > 0)
    )

    decreasing = bool(
        np.all(differences < 0)
    )

    if not (increasing or decreasing):
        raise ValueError(
            f"{axis_name} is not monotonic."
        )

    grid_index = np.arange(
        values.size,
        dtype=np.float64,
    )

    slope, intercept = np.polyfit(
        grid_index,
        values,
        deg=1,
    )

    fitted_values = (
        intercept
        + slope * grid_index
    )

    resolution = abs(
        float(slope)
    )

    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError(
            f"Could not calculate a valid resolution for {axis_name}."
        )

    residuals = (
        values
        - fitted_values
    )

    maximum_residual = float(
        np.max(
            np.abs(residuals)
        )
    )

    residual_fraction = (
        maximum_residual
        / resolution
    )

    # The diagnostics showed residual fractions of approximately:
    #
    # rlon: 0.076
    # rlat: 0.140
    #
    # These are consistent across all tested models and result from
    # coordinates rounded to two decimal places. A limit of 0.25 cells
    # permits normal rounding but still rejects genuinely distorted axes.
    maximum_allowed_fraction = 0.25

    if residual_fraction > maximum_allowed_fraction:
        raise ValueError(
            f"{axis_name} cannot be represented reliably as a regular "
            f"grid. Estimated resolution={resolution:.12f}; "
            f"maximum residual={maximum_residual:.12f}; "
            f"residual fraction={residual_fraction:.3%}."
        )

    return (
        fitted_values,
        resolution,
        maximum_residual,
    )


def _get_transform(
    data: xr.DataArray,
    *,
    y_dim: str,
    x_dim: str,
) -> tuple[object, bool]:
    """
    Create a GeoTIFF affine transform for a NARCliM rotated grid.

    The saved rlon and rlat coordinates may appear irregular because they
    are rounded to two decimal places. This function reconstructs their
    underlying regular axes using linear fitting before calculating the
    raster transform.

    Parameters
    ----------
    data
        Two-dimensional spatial DataArray, or a DataArray containing the
        requested spatial dimensions.

    y_dim
        Name of the vertical spatial dimension, normally ``rlat``.

    x_dim
        Name of the horizontal spatial dimension, normally ``rlon``.

    Returns
    -------
    transform
        Rasterio affine transform.

    flip_y
        True when the array must be flipped vertically before writing so
        raster rows are stored from north to south.
    """

    if x_dim not in data.coords:
        raise ValueError(
            f"Horizontal coordinate '{x_dim}' was not found."
        )

    if y_dim not in data.coords:
        raise ValueError(
            f"Vertical coordinate '{y_dim}' was not found."
        )

    x_values = np.asarray(
        data.coords[x_dim].values,
        dtype=np.float64,
    )

    y_values = np.asarray(
        data.coords[y_dim].values,
        dtype=np.float64,
    )

    (
        fitted_x,
        x_resolution,
        x_max_residual,
    ) = _fit_regular_axis(
        x_values,
        axis_name=x_dim,
    )

    (
        fitted_y,
        y_resolution,
        y_max_residual,
    ) = _fit_regular_axis(
        y_values,
        axis_name=y_dim,
    )

    # Cell coordinates represent pixel centres. The raster origin must
    # therefore be shifted half a cell to the outer northwest corner.
    west = float(
        np.min(fitted_x)
        - x_resolution / 2
    )

    north = float(
        np.max(fitted_y)
        + y_resolution / 2
    )

    transform = from_origin(
        west,
        north,
        x_resolution,
        y_resolution,
    )

    # GeoTIFF rows run north to south. NARCliM rlat normally increases
    # from south to north, so the array must be flipped vertically.
    flip_y = bool(
        fitted_y[0]
        < fitted_y[-1]
    )

    print(
        "[RASTER TRANSFORM] "
        f"{x_dim}_resolution={x_resolution:.10f}; "
        f"{y_dim}_resolution={y_resolution:.10f}; "
        f"{x_dim}_max_rounding_residual={x_max_residual:.10f}; "
        f"{y_dim}_max_rounding_residual={y_max_residual:.10f}; "
        f"flip_y={flip_y}",
        flush=True,
    )

    return transform, flip_y


def _write_geotiff(
    *,
    data: xr.DataArray,
    source_ds: xr.Dataset,
    output_path: Path,
    variable: str | None = None,
    scenario: str | None = None,
    time_horizon: str | None = None,
    ensemble_stat: str | None = None,
    nodata: float = -9999.0,
    overwrite: bool = False,
) -> Path:
    """
    Write a NARCliM climate layer as a correctly georeferenced EPSG:4326
    GeoTIFF and preserve climate-period metadata in the raster tags.

    NARCliM's native rlon/rlat grid is curvilinear in geographic space.
    The true 2-D lon/lat coordinates are therefore used by the existing
    rasterisation workflow before writing.

    Stored raster tags include, where available:

        variable
        scenario
        time_horizon
        time_horizon_label
        start_year
        end_year
        mid_year
        ensemble_stat

    This makes downstream modules such as spatial_summary.py self-contained:
    users do not need to define TIME_WINDOWS again.
    """

    if (
        output_path.exists()
        and not overwrite
    ):
        return output_path

    # ==================================================================
    # REGRID FROM NARCLIM CURVILINEAR GRID TO EPSG:4326
    # ==================================================================

    array, transform = (
        _rasterize_curvilinear_to_wgs84(
            data=data,
            source_ds=source_ds,
            nodata=nodata,
        )
    )

    array = np.asarray(
        array,
        dtype=np.float32,
    )

    array = np.where(
        np.isfinite(array),
        array,
        nodata,
    ).astype(
        np.float32
    )

    # ==================================================================
    # ROUND EXPORTED CLIMATE VALUES
    #
    # Calculations remain at full precision in memory.
    # Most variables use 1 decimal; StormDaysGT* uses 3 decimals so rare
    # event frequencies such as 0.033 are not collapsed to 0.0.
    # ==================================================================
    valid_output_mask = (
        array != nodata
    )

    array[
        valid_output_mask
    ] = np.round(
        array[
            valid_output_mask
        ],
        _output_decimals(
            variable
        ),
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_path = (
        output_path.with_suffix(
            ".tmp.tif"
        )
    )

    temporary_path.unlink(
        missing_ok=True
    )

    # ==================================================================
    # TEMPORAL METADATA
    # ==================================================================

    period_metadata = None

    if time_horizon is not None:
        period_metadata = (
            _get_time_period_metadata(
                source_ds=source_ds,
                window=time_horizon,
            )
        )

    # ==================================================================
    # WRITE A NORMAL GEOGRAPHIC RASTER
    # ==================================================================

    with rasterio.open(
        temporary_path,
        "w",
        driver="GTiff",
        height=array.shape[0],
        width=array.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=nodata,
        compress="deflate",
        predictor=3,
        tiled=True,
        BIGTIFF="IF_SAFER",
    ) as destination:

        destination.write(
            array,
            1,
        )

        destination.set_band_description(
            1,
            data.name or "climate_index",
        )

        tags: dict[str, str] = {
            "source_variable":
                str(
                    data.name
                    or "climate_index"
                ),

            "source_grid":
                "NARCliM native curvilinear rlat/rlon grid",

            "output_crs":
                "EPSG:4326",

            "regridding_method":
                (
                    "geographic source-cell polygon rasterisation; "
                    "all_touched=True"
                ),
        }

        if variable is not None:
            tags[
                "variable"
            ] = str(
                variable
            )

        if scenario is not None:
            tags[
                "scenario"
            ] = str(
                scenario
            )

        if time_horizon is not None:
            tags[
                "time_horizon"
            ] = str(
                time_horizon
            )

        if ensemble_stat is not None:
            tags[
                "ensemble_stat"
            ] = str(
                ensemble_stat
            )

        if period_metadata is not None:

            tags[
                "time_horizon_label"
            ] = str(
                period_metadata[
                    "time_horizon_label"
                ]
            )

            for field in (
                "start_year",
                "end_year",
                "mid_year",
            ):

                value = period_metadata[
                    field
                ]

                if value is not None:
                    tags[
                        field
                    ] = str(
                        int(
                            value
                        )
                    )

        destination.update_tags(
            **tags
        )

    temporary_path.replace(
        output_path
    )

    return output_path


def _find_subset_files(
    *,
    input_root: Path,
    variable: str,
    scenario: str,
    window: str,
    gcms: Iterable[str] | None,
    rcms: Iterable[str] | None,
    feature_ids: Iterable[str] | None,
) -> list[Path]:
    base = input_root / variable / scenario / window

    if not base.exists():
        return []

    files = sorted(base.rglob("*.nc"))

    if gcms is not None:
        allowed_gcms = set(gcms)
        files = [
            path
            for path in files
            if any(part in allowed_gcms for part in path.parts)
        ]

    if rcms is not None:
        aliases = {
            "NARCliM2-0-WRF412R3": {"NARCliM2-0-WRF412R3", "R3"},
            "NARCliM2-0-WRF412R5": {"NARCliM2-0-WRF412R5", "R5"},
            "R3": {"NARCliM2-0-WRF412R3", "R3"},
            "R5": {"NARCliM2-0-WRF412R5", "R5"},
        }

        allowed_parts: set[str] = set()
        for rcm in rcms:
            allowed_parts.update(aliases.get(rcm, {rcm}))

        files = [
            path
            for path in files
            if any(part in allowed_parts for part in path.parts)
        ]

    if feature_ids is not None:
        allowed_ids = set(feature_ids)
        files = [
            path
            for path in files
            if path.parent.name in allowed_ids
        ]

    return files


def _parse_model_metadata(
    *,
    path: Path,
    input_root: Path,
) -> tuple[str, str, str]:
    relative = path.relative_to(input_root)
    parts = relative.parts

    if len(parts) < 7:
        raise ValueError(f"Unexpected workflow path structure: {path}")

    return parts[3], parts[4], parts[5]


def _combine_source_files_for_model(
    *,
    paths: list[Path],
    variable: str,
) -> tuple[xr.Dataset, xr.DataArray]:
    """
    Combine all local NetCDF subset files belonging to one GCM x RCM
    ensemble member and calculate the temporal mean across the selected
    climate horizon.

    Examples
    --------
    TXge35:
        Mean annual number of hot days across the selected years.

    R20mm:
        Mean annual number of >=20 mm rainfall days.

    FFDIgt50:
        Mean annual number of FFDI >50 days.

    StormDaysGT89:
        Mean annual number of days with daily maximum near-surface wind
        speed >89 km/h.

    SPI12:
        Mean across the available monthly SPI12 values.

    Notes
    -----
    `reference_ds` retains the complete spatial variables and CRS metadata
    from the first source file.

    The full multi-file time period is stored in reference_ds.attrs rather
    than replacing its time coordinate. This avoids removing the original
    climate variable from the reference dataset.
    """

    if not paths:
        raise ValueError(
            "At least one source file is required."
        )

    # ------------------------------------------------------------------
    # IMPORTANT:
    #
    # decode_timedelta=False prevents xarray from interpreting variables
    # such as StormDaysGT89 with units='days' as timedelta64 data.
    # The values are climate-index counts, not time durations.
    # ------------------------------------------------------------------

    datasets = [
        xr.open_dataset(
            path,
            decode_times=True,
            decode_timedelta=False,
        )
        for path in paths
    ]

    try:

        # ==============================================================
        # EXTRACT CLIMATE VARIABLE FROM EACH SOURCE FILE
        # ==============================================================

        arrays = [
            ds[
                _choose_data_variable(
                    ds,
                    variable,
                )
            ]
            for ds in datasets
        ]

        # ==============================================================
        # COMBINE FILES THROUGH TIME
        # ==============================================================

        if len(arrays) == 1:

            combined = arrays[0]

        else:

            time_names = {
                _find_time_name(array)
                for array in arrays
            }

            if (
                len(time_names) != 1
                or None in time_names
            ):
                raise ValueError(
                    "Multiple source files do not share "
                    "one common time coordinate."
                )

            time_name = next(
                iter(time_names)
            )

            # ----------------------------------------------------------
            # IMPORTANT:
            #
            # `arrays` contains DataArray objects.
            #
            # Therefore DO NOT use:
            #
            #     data_vars="minimal"
            #
            # because `data_vars` is only valid when concatenating
            # Dataset objects.
            # ----------------------------------------------------------

            combined = xr.concat(
                arrays,
                dim=time_name,
                join="exact",
            )

            combined = combined.sortby(
                time_name
            )

        # ==============================================================
        # CALCULATE TEMPORAL MEAN
        # ==============================================================

        time_name = _find_time_name(
            combined
        )

        if time_name is None:

            model_mean = combined

        else:

            model_mean = combined.mean(
                dim=time_name,
                skipna=True,
                keep_attrs=True,
            )

        model_mean = (
            model_mean
            .squeeze(drop=True)
            .load()
        )

        model_mean.name = variable

        # ==============================================================
        # REFERENCE DATASET
        #
        # Keep the original first dataset intact.
        # Do NOT remove its climate variable.
        # ==============================================================

        reference_ds = datasets[0]

        # ==============================================================
        # STORE FULL TIME-HORIZON METADATA
        #
        # The first source file may only represent one year.
        # Therefore derive the full period from the combined time
        # coordinate and save it in global attributes.
        # ==============================================================

        if (
            time_name is not None
            and time_name in combined.coords
        ):

            combined_time = combined[
                time_name
            ]

            try:

                years = np.asarray(
                    combined_time.dt.year.values,
                    dtype=int,
                )

                years = years[
                    np.isfinite(years)
                ]

            except Exception:

                years = np.array(
                    [],
                    dtype=int,
                )

            if years.size:

                start_year = int(
                    np.min(years)
                )

                end_year = int(
                    np.max(years)
                )

                mid_year = int(
                    round(
                        (
                            start_year
                            + end_year
                        )
                        / 2
                    )
                )

                reference_ds.attrs[
                    "combined_start_year"
                ] = start_year

                reference_ds.attrs[
                    "combined_end_year"
                ] = end_year

                reference_ds.attrs[
                    "combined_mid_year"
                ] = mid_year

        return (
            reference_ds,
            model_mean,
        )

    except Exception:

        # If processing fails, make sure all opened datasets are closed.
        for ds in datasets:
            ds.close()

        raise

def _calculate_ensemble_statistics(
    model_layers: list[ModelLayer],
) -> dict[str, xr.DataArray]:
    """
    Calculate pixel-wise ensemble statistics across GCM × RCM layers.
    """

    if not model_layers:
        raise ValueError(
            "No model layers were supplied."
        )

    reference = model_layers[0].data

    aligned_layers = []
    ensemble_labels = []

    for layer in model_layers:
        reference_aligned, candidate_aligned = xr.align(
            reference,
            layer.data,
            join="exact",
            copy=False,
        )

        # This should normally remain unchanged, but using the aligned
        # reference makes the intent explicit.
        reference = reference_aligned

        aligned_layers.append(
            candidate_aligned
        )

        ensemble_labels.append(
            f"{layer.gcm}_{layer.rcm}"
        )

    ensemble = xr.concat(
        aligned_layers,
        dim=pd.Index(
            ensemble_labels,
            name="ensemble_member",
        ),
    )

    statistics = {
        "mean": ensemble.mean(
            dim="ensemble_member",
            skipna=True,
            keep_attrs=True,
        ),
        "min": ensemble.min(
            dim="ensemble_member",
            skipna=True,
            keep_attrs=True,
        ),
        "max": ensemble.max(
            dim="ensemble_member",
            skipna=True,
            keep_attrs=True,
        ),
    }

    for statistic_name, statistic_data in statistics.items():
        statistic_data.attrs.update(
            {
                "ensemble_statistic": statistic_name,
                "ensemble_member_count": len(
                    model_layers
                ),
                "ensemble_members": ", ".join(
                    ensemble_labels
                ),
            }
        )

    return statistics



def _get_time_period_metadata(
    *,
    source_ds: xr.Dataset,
    window: str,
) -> dict[str, object]:
    """
    Read the ACTUAL temporal period represented by a local NetCDF subset.

    The uncertainty workflow should not maintain a second TIME_WINDOWS
    dictionary. The local NetCDF produced by workflow.py already contains the
    actual time coordinate for the requested horizon, so this function derives
    the period directly from that coordinate.

    Returns
    -------
    dict
        Example:
        {
            "time_horizon": "mid_term",
            "start_year": 2041,
            "end_year": 2060,
            "mid_year": 2050,
            "time_horizon_label": "mid_term (2050)",
        }

    Notes
    -----
    For an even-length interval there are two central calendar years. The
    lower integer midpoint is used, matching the package's earlier convention:
    2021-2040 -> 2030.
    """
    # ==================================================================
    # PREFER FULL MULTI-FILE PERIOD METADATA
    #
    # _combine_source_files_for_model() stores these attributes after
    # combining all annual/monthly source files belonging to a climate
    # horizon.
    #
    # This is important when each local NetCDF contains only one year.
    # ==================================================================

    combined_start_year = (
        source_ds.attrs.get(
            "combined_start_year"
        )
    )

    combined_end_year = (
        source_ds.attrs.get(
            "combined_end_year"
        )
    )

    combined_mid_year = (
        source_ds.attrs.get(
            "combined_mid_year"
        )
    )

    if (
        combined_start_year is not None
        and combined_end_year is not None
    ):

        start_year = int(
            combined_start_year
        )

        end_year = int(
            combined_end_year
        )

        if combined_mid_year is None:

            mid_year = int(
                round(
                    (
                        start_year
                        + end_year
                    )
                    / 2
                )
            )

        else:

            mid_year = int(
                combined_mid_year
            )

        return {
            "time_horizon":
                window,

            "time_horizon_label":
                (
                    f"{window} "
                    f"({start_year}-{end_year}; "
                    f"mid={mid_year})"
                ),

            "start_year":
                start_year,

            "end_year":
                end_year,

            "mid_year":
                mid_year,
        }
    metadata: dict[str, object] = {
        "time_horizon": str(window),
        "start_year": None,
        "end_year": None,
        "mid_year": None,
        "time_horizon_label": str(window),
    }

    time_name = _find_time_name(
        source_ds
    )

    if time_name is None:
        return metadata

    try:
        years = np.asarray(
            source_ds[
                time_name
            ].dt.year.values,
            dtype=np.int64,
        )
    except Exception:
        return metadata

    years = years[
        np.isfinite(
            years
        )
    ]

    if years.size == 0:
        return metadata

    start_year = int(
        np.nanmin(
            years
        )
    )

    end_year = int(
        np.nanmax(
            years
        )
    )

    mid_year = (
        start_year
        + end_year
    ) // 2

    metadata.update(
        {
            "start_year":
                start_year,

            "end_year":
                end_year,

            "mid_year":
                mid_year,

            "time_horizon_label":
                f"{window} ({mid_year})",
        }
    )

    return metadata


def _format_time_horizon(
    *,
    source_ds: xr.Dataset,
    window: str,
) -> str:
    """
    Return a readable horizon label using the actual NetCDF time coordinate.

    Example:
        mid_term -> mid_term (2050)
    """

    return str(
        _get_time_period_metadata(
            source_ds=source_ds,
            window=window,
        )[
            "time_horizon_label"
        ]
    )


def _raster_pixels_to_points(
    *,
    raster_path: Path,
    feature_id: str,
    variable: str,
    scenario: str,
    time_horizon: str,
    start_year: int | None,
    end_year: int | None,
    mid_year: int | None,
    gcm: str,
    rcm: str,
) -> gpd.GeoDataFrame:
    """
    Convert every valid raster pixel to a point at the pixel centre.

    This layer represents the temporal-average value of each individual
    GCM x RCM member. It is useful for QA and for inspecting the model
    values that contribute to the ensemble statistics.
    """

    with rasterio.open(raster_path) as src:

        if src.crs is None:
            raise ValueError(
                f"Raster has no CRS: {raster_path}"
            )

        raster = src.read(
            1,
            masked=True,
        )

        valid_mask = ~np.ma.getmaskarray(
            raster
        )

        rows, cols = np.where(
            valid_mask
        )

        if rows.size == 0:
            return gpd.GeoDataFrame(
                columns=[
                    "pixel_id",
                    "feature_id",
                    "variable",
                    "scenario",
                    "time_horizon",
                    "start_year",
                    "end_year",
                    "mid_year",
                    "gcm",
                    "rcm",
                    "value",
                    "longitude",
                    "latitude",
                    "raster",
                    "geometry",
                ],
                geometry="geometry",
                crs=src.crs,
            )

        xs, ys = rasterio.transform.xy(
            src.transform,
            rows,
            cols,
            offset="center",
        )

        xs = np.asarray(
            xs,
            dtype=np.float64,
        )

        ys = np.asarray(
            ys,
            dtype=np.float64,
        )

        values = np.asarray(
            raster.data[
                rows,
                cols,
            ],
            dtype=np.float64,
        )

        # Keep point attributes consistent with raster precision.
        values = np.round(
            values,
            _output_decimals(
                variable
            ),
        )

        # Coordinates are used in the ID because all output rasters are
        # snapped to the package's fixed EPSG:4326 analysis grid.
        pixel_ids = [
            f"{x:.6f}_{y:.6f}"
            for x, y in zip(
                xs,
                ys,
            )
        ]

        geometry = [
            Point(
                float(x),
                float(y),
            )
            for x, y in zip(
                xs,
                ys,
            )
        ]

        return gpd.GeoDataFrame(
            {
                "pixel_id":
                    pixel_ids,

                "feature_id":
                    str(feature_id),

                "variable":
                    variable,

                "scenario":
                    scenario,

                "time_horizon":
                    time_horizon,

                "start_year":
                    start_year,

                "end_year":
                    end_year,

                "mid_year":
                    mid_year,

                "gcm":
                    gcm,

                "rcm":
                    rcm,

                "value":
                    values,

                "longitude":
                    xs,

                "latitude":
                    ys,

                "raster":
                    raster_path.name,

                "geometry":
                    geometry,
            },
            geometry="geometry",
            crs=src.crs,
        )


def _ensemble_rasters_to_points(
    *,
    raster_paths: dict[str, Path],
    feature_id: str,
    variable: str,
    scenario: str,
    time_horizon: str,
    start_year: int | None,
    end_year: int | None,
    mid_year: int | None,
) -> gpd.GeoDataFrame:
    """
    Create one point per valid climate pixel containing the pixel-wise
    ensemble mean, minimum and maximum.

    The output contains NO GCM or RCM fields because the purpose of this
    layer is to represent the uncertainty envelope across the full model
    ensemble, rather than individual model values.
    """

    required_statistics = {
        "mean",
        "min",
        "max",
    }

    missing = (
        required_statistics
        - set(
            raster_paths
        )
    )

    if missing:
        raise ValueError(
            "Missing ensemble raster(s): "
            f"{sorted(missing)}"
        )

    arrays = {}
    reference_profile = None
    reference_transform = None
    reference_crs = None
    reference_shape = None

    for statistic in (
        "mean",
        "min",
        "max",
    ):

        raster_path = raster_paths[
            statistic
        ]

        with rasterio.open(
            raster_path
        ) as src:

            if src.crs is None:
                raise ValueError(
                    f"Raster has no CRS: {raster_path}"
                )

            current_shape = (
                src.height,
                src.width,
            )

            if reference_profile is None:

                reference_profile = (
                    src.profile.copy()
                )

                reference_transform = (
                    src.transform
                )

                reference_crs = (
                    src.crs
                )

                reference_shape = (
                    current_shape
                )

            else:

                if current_shape != reference_shape:
                    raise ValueError(
                        "Ensemble rasters do not share the same shape. "
                        f"Expected {reference_shape}; "
                        f"{statistic} has {current_shape}."
                    )

                if not src.transform.almost_equals(
                    reference_transform
                ):
                    raise ValueError(
                        "Ensemble rasters do not share the same "
                        "affine transform."
                    )

                if src.crs != reference_crs:
                    raise ValueError(
                        "Ensemble rasters do not share the same CRS."
                    )

            arrays[
                statistic
            ] = src.read(
                1,
                masked=True,
            )

    # A point is retained where at least one ensemble statistic is valid.
    valid_mask = np.zeros(
        reference_shape,
        dtype=bool,
    )

    for array in arrays.values():
        valid_mask |= (
            ~np.ma.getmaskarray(
                array
            )
        )

    rows, cols = np.where(
        valid_mask
    )

    if rows.size == 0:
        return gpd.GeoDataFrame(
            columns=[
                "feature_id",
                "variable",
                "scenario",
                "time_horizon",
                "start_year",
                "end_year",
                "mid_year",
                "mean",
                "min",
                "max",
                "longitude",
                "latitude",
                "geometry",
            ],
            geometry="geometry",
            crs=reference_crs,
        )

    xs, ys = rasterio.transform.xy(
        reference_transform,
        rows,
        cols,
        offset="center",
    )

    xs = np.asarray(
        xs,
        dtype=np.float64,
    )

    ys = np.asarray(
        ys,
        dtype=np.float64,
    )

    def _values_for(
        statistic: str,
    ) -> np.ndarray:
        array = arrays[
            statistic
        ]

        values = np.asarray(
            array.data[
                rows,
                cols,
            ],
            dtype=np.float64,
        )

        masked = np.ma.getmaskarray(
            array
        )[
            rows,
            cols,
        ]

        values[
            masked
        ] = np.nan

        # Keep ensemble point attributes consistent with raster precision.
        values = np.round(
            values,
            _output_decimals(
                variable
            ),
        )

        return values


    geometry = [
        Point(
            float(x),
            float(y),
        )
        for x, y in zip(
            xs,
            ys,
        )
    ]

    return gpd.GeoDataFrame(
        {

            "feature_id":
                str(feature_id),

            "variable":
                variable,

            "scenario":
                scenario,

            "time_horizon":
                time_horizon,

            "start_year":
                start_year,

            "end_year":
                end_year,

            "mid_year":
                mid_year,

            # These three fields are the pixel-wise model-ensemble
            # uncertainty statistics requested for the final summary.
            "mean":
                _values_for(
                    "mean"
                ),

            "min":
                _values_for(
                    "min"
                ),

            "max":
                _values_for(
                    "max"
                ),


            "longitude":
                xs,

            "latitude":
                ys,

            "geometry":
                geometry,
        },
        geometry="geometry",
        crs=reference_crs,
    )


def _summarise_raster_for_boundaries(
    *,
    raster_path: Path,
    boundaries: gpd.GeoDataFrame,
    variable: str,
    scenario: str,
    window: str,
    statistic: str,
    ensemble_members: int,
    id_field: str | None,
) -> gpd.GeoDataFrame:
    """
    Calculate spatial statistics from an ensemble GeoTIFF for each
    supplied boundary.

    The uncertainty GeoTIFFs are written in EPSG:4326, so the supplied
    boundaries are transformed into the raster CRS before masking.
    """

    if boundaries.empty:
        raise ValueError(
            "No boundary features were supplied."
        )

    if boundaries.crs is None:
        raise ValueError(
            "The boundary layer has no CRS."
        )

    rows = []

    with rasterio.open(
        raster_path
    ) as src:

        if src.crs is None:
            raise ValueError(
                f"Raster has no CRS: {raster_path}"
            )

        boundary_data = (
            boundaries.to_crs(
                src.crs
            )
        )

        raster = src.read(
            1,
            masked=True,
        )

        for position, (
            (_, original_row),
            (_, raster_row),
        ) in enumerate(
            zip(
                boundaries.iterrows(),
                boundary_data.iterrows(),
            ),
            start=1,
        ):

            geometry = (
                raster_row.geometry
            )

            if (
                geometry is None
                or geometry.is_empty
            ):

                values = np.array(
                    [],
                    dtype=np.float64,
                )

            else:

                inside_boundary = geometry_mask(
                    [
                        geometry
                    ],

                    out_shape=(
                        src.height,
                        src.width,
                    ),

                    transform=(
                        src.transform
                    ),

                    invert=True,

                    # Include any raster cell touched by the boundary.
                    all_touched=True,
                )

                combined_mask = (
                    np.ma.getmaskarray(
                        raster
                    )
                    |
                    ~inside_boundary
                )

                values = np.ma.array(
                    raster,
                    mask=combined_mask,
                ).compressed()

            # ==========================================================
            # STATISTICS
            # ==========================================================

            if values.size:

                spatial_mean = round(
                    float(
                        np.nanmean(values)
                    ),
                    1,
                )

                spatial_min = round(
                    float(
                        np.nanmin(values)
                    ),
                    1,
                )

                spatial_max = round(
                    float(
                        np.nanmax(values)
                    ),
                    1,
                )

                valid_cell_count = int(
                    values.size
                )

            else:

                spatial_mean = np.nan
                spatial_min = np.nan
                spatial_max = np.nan

                valid_cell_count = 0

            # ==========================================================
            # FEATURE ID
            # ==========================================================

            if (
                id_field
                and id_field
                in original_row.index
            ):

                feature_id = (
                    original_row[
                        id_field
                    ]
                )

            else:

                feature_id = (
                    f"ID{position}"
                )

            rows.append(
                {
                    "feature_id":
                        str(
                            feature_id
                        ),

                    "variable":
                        variable,

                    "scenario":
                        scenario,

                    "time_horizon":
                        window,

                    "ensemble_stat":
                        statistic,

                    "mean":
                        spatial_mean,

                    "min":
                        spatial_min,

                    "max":
                        spatial_max,

                    "valid_cells":
                        valid_cell_count,

                    "n_models":
                        ensemble_members,

                    "raster":
                        raster_path.name,

                    # Keep original geometry.
                    "geometry":
                        original_row.geometry,
                }
            )

    return gpd.GeoDataFrame(
        rows,
        geometry="geometry",
        crs=boundaries.crs,
    )


# =============================================================================
# VARIABLE AND SCENARIO AVAILABILITY
# =============================================================================

def _available_scenarios_for_variable(
    variable: str,
) -> tuple[str, ...]:
    """
    Return the scenarios configured for a climate variable.

    Older VariableSpec definitions that do not contain
    ``available_scenarios`` remain compatible and default to all standard
    NARCliM scenarios.
    """

    if variable in VARIABLES:
        specification = VARIABLES[
            variable
        ]

    elif _is_storm_derived_variable(
        variable
    ):
        # Dynamically named storm outputs are generated by workflow.py
        # from sfcWindmax using a user-defined threshold.
        return (
            "historical",
            "ssp126",
            "ssp245",
            "ssp370",
        )

    else:
        raise KeyError(
            f"Variable '{variable}' is not configured. "
            f"Available variables: {sorted(VARIABLES)}"
        )

    return tuple(
        getattr(
            specification,
            "available_scenarios",
            (
                "historical",
                "ssp126",
                "ssp245",
                "ssp370",
            ),
        )
    )


def _scenario_window_is_valid(
    *,
    scenario: str,
    window: str,
) -> tuple[bool, str]:
    """
    Validate the standard scenario/time-window combinations.

    This package convention pairs:

    - ``historical`` with ``baseline``;
    - SSP scenarios with future windows such as ``short_term``,
      ``mid_term`` and ``long_term``.

    Other custom future-window names are allowed as long as they are not
    named ``baseline``.
    """

    if scenario == "historical" and window != "baseline":
        return (
            False,
            "Historical data are only paired with the baseline window.",
        )

    if scenario != "historical" and window == "baseline":
        return (
            False,
            "The baseline window is paired with historical data, "
            "not projection scenarios.",
        )

    return True, ""

def run_uncertainty_workflow(
    *,
    input_root: str | Path,
    boundary_path: str | Path,
    variables: list[str],
    scenarios: list[str],
    time_windows: list[str],
    gcms: list[str] | None = None,
    rcms: list[str] | None = None,
    feature_ids: list[str] | None = None,
    boundary_id_field: str | None = None,
    output_folder_name: str = "Climate_Indices",
    overwrite: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:

    # ==================================================================
    # PROJECT, DATA, AND OUTPUT FOLDERS
    # ==================================================================

    # Main project/output root supplied by the user.
    #
    # Example:
    #   C:\NARCliM_Outputs
    input_root = Path(input_root)

    # Downloaded/subset NARCliM NetCDF data created by workflow.py.
    #
    # Example:
    #   C:\NARCliM_Outputs\NARCliM Data
    data_root = ( input_root / "NARCliM Data")
    # Outputs created by uncertainty.py.
    #
    # Example:
    #   C:\NARCliM_Outputs\Climate_Indices
    output_root = (input_root / output_folder_name)
    # ------------------------------------------------------------------
    # Uncertainty output folders
    # ------------------------------------------------------------------
    temporal_average_folder = (output_root / "Temporal_Average")
    ensemble_folder = (output_root / "Ensemble_Stats")
    hazards_folder = (output_root / "Hazards")
    temporal_average_folder.mkdir(parents=True,exist_ok=True,)
    ensemble_folder.mkdir(parents=True,exist_ok=True,)
    hazards_folder.mkdir(parents=True,exist_ok=True,)
    # ------------------------------------------------------------------
    # Validate the NARCliM data folder
    # ------------------------------------------------------------------
    if not data_root.exists():
        raise FileNotFoundError(
            "NARCliM data folder was not found:\n"
            f"{data_root}\n\n"
            "Run workflow.py first or check input_root.")
    # ------------------------------------------------------------------
    # Read study boundary
    # ------------------------------------------------------------------
    boundaries = gpd.read_file(boundary_path)
    if boundaries.empty:
        raise ValueError(f"No boundary features found: {boundary_path}")
    if boundaries.crs is None:
        raise ValueError("Boundary file has no CRS.")
    # ------------------------------------------------------------------
    # Containers used throughout the workflow
    # ------------------------------------------------------------------
    records = []
    gpkg_model_points = []
    gpkg_ensemble_points = []

    _log("=" * 100,verbose=verbose,)
    _log("NARCliM CLIMATE-MODEL UNCERTAINTY WORKFLOW",verbose=verbose,)
    _log("=" * 100,verbose=verbose,)
    _log(f"[CODE VERSION] {UNCERTAINTY_CODE_VERSION}",verbose=verbose,)
    _log(f"[NARCliM DATA] {data_root}",verbose=verbose,)
    _log(f"[OUTPUT ROOT] {output_root}",verbose=verbose,)
  
    for variable in variables:
        available_scenarios = _available_scenarios_for_variable(variable)

        for scenario in scenarios:
            if scenario not in available_scenarios:
                message = (
                    f"{variable} is not available for scenario '{scenario}'."
                )
                _log(
                    f"[DATA NOT AVAILABLE] {message}",
                    verbose=verbose,
                )
                records.append({
                    "status": "data_not_available",
                    "variable": variable,
                    "scenario": scenario,
                    "time_window": None,
                    "gcm": None,
                    "rcm": None,
                    "feature_id": None,
                    "output": None,
                    "message": message,
                })
                continue

            for window in time_windows:
                _log(
                    f"[GROUP] {variable} | {scenario} | {window}",
                    verbose=verbose,
                )

                valid, reason = _scenario_window_is_valid(
                    scenario=scenario,
                    window=window,
                )

                if not valid:
                    _log(f"[SKIP] {reason}", verbose=verbose)
                    records.append({
                        "status": "time_window_not_available",
                        "variable": variable,
                        "scenario": scenario,
                        "time_window": window,
                        "gcm": None,
                        "rcm": None,
                        "feature_id": None,
                        "output": None,
                        "message": reason,
                    })
                    continue

                source_files = _find_subset_files(
                    input_root=data_root,
                    variable=variable,
                    scenario=scenario,
                    window=window,
                    gcms=gcms,
                    rcms=rcms,
                    feature_ids=feature_ids,
                )

                if not source_files:
                    records.append({
                        "status": "no_input_files",
                        "variable": variable,
                        "scenario": scenario,
                        "time_window": window,
                        "gcm": None,
                        "rcm": None,
                        "feature_id": None,
                        "output": None,
                        "message": "No input NetCDF subset files found.",
                    })
                    _log("[DATA NOT AVAILABLE] No input NetCDF files found.", verbose=verbose)
                    continue

                groups = {}

                for path in source_files:
                    gcm, rcm, feature_id = _parse_model_metadata(
                        path=path,
                        input_root=data_root,
                    )
                    groups.setdefault(
                        (gcm, rcm, feature_id),
                        [],
                    ).append(path)

                feature_layers = {}
                reference_datasets = {}

                for (gcm, rcm, feature_id), paths in sorted(groups.items()):
                    try:
                        source_ds, model_mean = _combine_source_files_for_model(
                            paths=sorted(paths),
                            variable=variable,
                        )

                        reference_datasets.setdefault(
                            feature_id,
                            source_ds,
                        )

                        output_name = (
                            f"{_safe_name(variable)}_"
                            f"{_safe_name(scenario)}_"
                            f"{_safe_name(window)}_"
                            f"{_safe_name(gcm)}_"
                            f"{_safe_name(rcm)}_"
                            f"{_safe_name(feature_id)}_"
                            f"temporal_mean.tif"
                        )

                        output_path = temporal_average_folder / output_name

                        period_metadata = _get_time_period_metadata(
                            source_ds=source_ds,
                            window=window,
                        )

                        _write_geotiff(
                            data=model_mean,
                            source_ds=source_ds,
                            output_path=output_path,
                            variable=variable,
                            scenario=scenario,
                            time_horizon=window,
                            overwrite=overwrite,
                        )

                        gpkg_model_points.append(
                            _raster_pixels_to_points(
                                raster_path=output_path,
                                feature_id=feature_id,
                                variable=variable,
                                scenario=scenario,
                                time_horizon=str(
                                    period_metadata[
                                        "time_horizon_label"
                                    ]
                                ),
                                start_year=period_metadata[
                                    "start_year"
                                ],
                                end_year=period_metadata[
                                    "end_year"
                                ],
                                mid_year=period_metadata[
                                    "mid_year"
                                ],
                                gcm=gcm,
                                rcm=rcm,
                            )
                        )

                        feature_layers.setdefault(
                            feature_id,
                            [],
                        ).append(
                            ModelLayer(
                                variable=variable,
                                scenario=scenario,
                                window=window,
                                gcm=gcm,
                                rcm=rcm,
                                feature_id=feature_id,
                                source_path=paths[0],
                                data=model_mean,
                            )
                        )

                        records.append({
                            "status": "model_layer_saved",
                            "variable": variable,
                            "scenario": scenario,
                            "time_window": window,
                            "gcm": gcm,
                            "rcm": rcm,
                            "feature_id": feature_id,
                            "output": str(output_path),
                            "message": "",
                        })

                        _log(f"[OK] {output_path.name}", verbose=verbose)

                    except Exception as exc:
                        traceback.print_exc()
                        message = f"{type(exc).__name__}: {exc}"

                        records.append({
                            "status": "model_layer_error",
                            "variable": variable,
                            "scenario": scenario,
                            "time_window": window,
                            "gcm": gcm,
                            "rcm": rcm,
                            "feature_id": feature_id,
                            "output": None,
                            "message": message,
                        })

                        _log(f"[ERROR] {message}", verbose=verbose)

                for feature_id, model_layers in feature_layers.items():
                    if not model_layers:
                        continue

                    try:
                        stats = _calculate_ensemble_statistics(
                            model_layers
                        )

                        reference_ds = reference_datasets[
                            feature_id
                        ]

                        # Recover the authoritative NARCliM rotated-pole CRS from the
                        # source NetCDF. This is also used when summarising the GeoTIFF.
                        reference_variable_name = _choose_data_variable(
                            reference_ds,
                            variable,
                        )

                        reference_crs = _get_crs(
                            reference_ds,
                            reference_ds[reference_variable_name],
                        )

                        #_log(f"[ENSEMBLE CRS] {reference_crs.to_string()}",verbose=verbose,)

                        # --------------------------------------------------------------
                        # Determine the actual period represented by this ensemble.
                        # --------------------------------------------------------------

                        period_metadata = _get_time_period_metadata(
                            source_ds=reference_ds,
                            window=window,
                        )

                        start_year = period_metadata[
                            "start_year"
                        ]

                        end_year = period_metadata[
                            "end_year"
                        ]

                        mid_year = period_metadata[
                            "mid_year"
                        ]

                        if (
                            start_year is None
                            or end_year is None
                        ):
                            raise ValueError(
                                "Could not determine start/end year for "
                                f"{variable} | {scenario} | {window}."
                            )

                        # --------------------------------------------------------------
                        # Write the three ensemble-statistic rasters using clean names.
                        #
                        # Examples:
                        #   TXge35_Min_historical_1985_2014.tif
                        #   TXge35_Mean_ssp245_2040_2060.tif
                        #   TXge35_Max_ssp370_2080_2100.tif
                        # --------------------------------------------------------------

                        ensemble_raster_paths = {}

                        for statistic, statistic_data in stats.items():

                            statistic_data.name = variable

                            statistic_label = (
                                statistic.capitalize()
                            )

                            # Keep the requested clean filename for the normal
                            # single-feature workflow. Additional features retain
                            # their ID so outputs cannot overwrite one another.
                            feature_suffix = (
                                ""
                                if str(feature_id) == "ID1"
                                else f"_{_safe_name(feature_id)}"
                            )

                            output_name = (
                                f"{_safe_name(variable)}_"
                                f"{statistic_label}_"
                                f"{_safe_name(scenario)}_"
                                f"{int(start_year)}_"
                                f"{int(end_year)}"
                                f"{feature_suffix}.tif"
                            )

                            output_path = (
                                ensemble_folder
                                / output_name
                            )

                            _write_geotiff(
                                data=statistic_data,
                                source_ds=reference_ds,
                                output_path=output_path,
                                variable=variable,
                                scenario=scenario,
                                time_horizon=window,
                                ensemble_stat=statistic,
                                overwrite=overwrite,
                            )

                            records.append({
                                "status": "ensemble_raster_saved",
                                "variable": variable,
                                "scenario": scenario,
                                "time_window": window,
                                "gcm": None,
                                "rcm": None,
                                "feature_id": feature_id,
                                "output": str(output_path),
                                "message": "",
                            })

                            ensemble_raster_paths[
                                statistic
                            ] = output_path

                            _log(
                                f"[OK] {output_path.name}",
                                verbose=verbose,
                            )

                        gpkg_ensemble_points.append(
                            _ensemble_rasters_to_points(
                                raster_paths=ensemble_raster_paths,
                                feature_id=feature_id,
                                variable=variable,
                                scenario=scenario,
                                time_horizon=str(
                                    period_metadata[
                                        "time_horizon_label"
                                    ]
                                ),
                                start_year=period_metadata[
                                    "start_year"
                                ],
                                end_year=period_metadata[
                                    "end_year"
                                ],
                                mid_year=period_metadata[
                                    "mid_year"
                                ],
                            )
                        )

                    except Exception as exc:
                        traceback.print_exc()
                        message = f"{type(exc).__name__}: {exc}"

                        records.append({
                            "status": "ensemble_error",
                            "variable": variable,
                            "scenario": scenario,
                            "time_window": window,
                            "gcm": None,
                            "rcm": None,
                            "feature_id": feature_id,
                            "output": None,
                            "message": message,
                        })

                        _log(f"[ENSEMBLE ERROR] {message}", verbose=verbose)

                for ds in reference_datasets.values():
                    ds.close()

    # ==================================================================
    # HAZARD GEOPACKAGES
    #
    # The final pixel-wise ensemble hazard products are written to:
    #
    #   Climate_Indices/
    #       Hazards/
    #           Ensemble_summary.gpkg
    #           <VARIABLE>.gpkg
    #
    # Ensemble_summary.gpkg contains all requested variables in one
    # point layer named ``pixel_ensemble_summary``.
    #
    # Each variable GeoPackage contains only that variable's points and
    # uses the variable name itself as the layer name.
    #
    # The values are:
    #
    #   1. temporal mean across all years/months in the selected horizon
    #      for each individual GCM x RCM member;
    #   2. pixel-wise ensemble mean/min/max across those model-member
    #      temporal means.
    # ==================================================================

    has_ensemble_points = any(
        not gdf.empty
        for gdf in gpkg_ensemble_points
    )

    if has_ensemble_points:

        valid_ensemble_frames = [
            gdf
            for gdf in gpkg_ensemble_points
            if not gdf.empty
        ]

        ensemble_points = pd.concat(
            valid_ensemble_frames,
            ignore_index=True,
        )

        ensemble_points = gpd.GeoDataFrame(
            ensemble_points,
            geometry="geometry",
            crs=valid_ensemble_frames[0].crs,
        )

        # --------------------------------------------------------------
        # Final attribute precision:
        #   most variables -> 1 decimal
        #   StormDaysGT*    -> 3 decimals
        # --------------------------------------------------------------
        for climate_field in (
            "mean",
            "min",
            "max",
            "value",
        ):
            if climate_field not in ensemble_points.columns:
                continue

            ensemble_points[
                climate_field
            ] = pd.to_numeric(
                ensemble_points[
                    climate_field
                ],
                errors="coerce",
            )

            for variable_name in (
                ensemble_points[
                    "variable"
                ]
                .dropna()
                .unique()
            ):
                variable_mask = (
                    ensemble_points[
                        "variable"
                    ]
                    == variable_name
                )

                ensemble_points.loc[
                    variable_mask,
                    climate_field,
                ] = (
                    ensemble_points.loc[
                        variable_mask,
                        climate_field,
                    ]
                    .round(
                        _output_decimals(
                            variable_name
                        )
                    )
                )

        # --------------------------------------------------------------
        # 1. COMBINED GPKG: all variables together.
        # --------------------------------------------------------------

        ensemble_summary_path = (
            hazards_folder
            / "Ensemble_summary_all_variables.gpkg"
        )

        if ensemble_summary_path.exists():
            ensemble_summary_path.unlink()

        ensemble_points.to_file(
            ensemble_summary_path,
            layer="ensemble_summary_all_variables",
            driver="GPKG",
            index=False,
        )

        _log(
            "[HAZARD GPKG] Ensemble_summary_all_variables.gpkg | "
            f"{len(ensemble_points)} point record(s)",
            verbose=verbose,
        )

        # --------------------------------------------------------------
        # 2. ONE GPKG PER VARIABLE.
        #
        # Example:
        #   Hazards/Ensemble_summary_TXge35.gpkg
        #   Hazards/Ensemble_summary_R20mm.gpkg
        #   Hazards/Ensemble_summary_R99p.gpkg
        #   Hazards/Ensemble_summary_FFDIgt50.gpkg
        #   Hazards/Ensemble_summary_SPI12.gpkg
        # --------------------------------------------------------------

        for variable_name in sorted(
            ensemble_points[
                "variable"
            ].dropna().unique()
        ):

            variable_points = (
                ensemble_points[
                    ensemble_points[
                        "variable"
                    ]
                    == variable_name
                ]
                .copy()
            )

            if variable_points.empty:
                continue

            safe_variable_name = (
                _safe_name(
                    variable_name
                )
            )

            variable_file_name = (
                f"Ensemble_summary_{safe_variable_name}.gpkg"
            )

            variable_gpkg_path = (
                hazards_folder
                / variable_file_name
            )

            if variable_gpkg_path.exists():
                variable_gpkg_path.unlink()

            variable_layer_name = (
                f"ensemble_summary_{safe_variable_name}"
            )

            variable_points.to_file(
                variable_gpkg_path,
                layer=variable_layer_name,
                driver="GPKG",
                index=False,
            )

            _log(
                "[HAZARD GPKG] "
                f"{variable_file_name} | "
                f"{len(variable_points)} point record(s)",
                verbose=verbose,
            )

        _log(
            f"\n[HAZARDS FOLDER] {hazards_folder}",
            verbose=verbose,
        )

    # --------------------------------------------------------------
    # UNCERTAINTY MANIFEST
    #
    # Always create the manifest, including runs where no ensemble
    # point GeoPackage was produced.
    # --------------------------------------------------------------

    manifest = pd.DataFrame(
        records
    )

    manifest_path = (
        output_root
        / "uncertainty_manifest.csv"
    )

    manifest.to_csv(
        manifest_path,
        index=False,
    )

    _log(
        f"\n[MANIFEST] {manifest_path}",
        verbose=verbose,
    )

    if not manifest.empty:
        _log(
            "\n[STATUS SUMMARY]\n"
            + manifest["status"].value_counts(dropna=False).to_string(),
            verbose=verbose,
        )

    return manifest