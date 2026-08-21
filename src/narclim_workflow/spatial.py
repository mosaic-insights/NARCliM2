from dataclasses import dataclass
from pathlib import Path
import re

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point, Polygon
from shapely.prepared import prep


LAT_CANDIDATES = ("lat", "latitude", "nav_lat", "rlat")
LON_CANDIDATES = ("lon", "longitude", "nav_lon", "rlon")
TIME_CANDIDATES = ("time",)


@dataclass(frozen=True)
class SpatialFeature:
    """One vector feature supplied to the NARCliM workflow."""

    feature_id: str
    geometry: object


def safe_name(value: object) -> str:
    """Convert a feature ID into a filesystem-safe string."""

    text = str(value).strip()
    text = re.sub(r'[<>:"/\\|?*]+', "_", text)
    text = re.sub(r"\s+", "_", text)

    return text.strip("._") or "feature"


def read_features(
    vector_path: str | Path,
    *,
    id_field: str | None = None,
) -> list[SpatialFeature]:
    """
    Read point/polygon features and convert them to EPSG:4326.

    When ``id_field`` is not supplied, features are named ID1, ID2, ...
    """

    gdf = gpd.read_file(vector_path)

    if gdf.empty:
        raise ValueError(
            f"No features found: {vector_path}"
        )

    if gdf.crs is None:
        raise ValueError(
            "The input vector file has no CRS."
        )

    # NCI's 2-D lon/lat coordinates are geographic longitude/latitude.
    gdf = gdf.to_crs("EPSG:4326")

    features: list[SpatialFeature] = []

    for position, (_, row) in enumerate(
        gdf.iterrows(),
        start=1,
    ):
        geometry = row.geometry

        if geometry is None or geometry.is_empty:
            continue

        if id_field is not None:

            if id_field not in gdf.columns:
                raise KeyError(
                    f"id_field '{id_field}' was not found. "
                    f"Available fields: {list(gdf.columns)}"
                )

            raw_id = row[id_field]

            if raw_id is None or str(raw_id).strip() == "":
                raw_id = f"ID{position}"

        else:
            raw_id = f"ID{position}"

        features.append(
            SpatialFeature(
                feature_id=safe_name(raw_id),
                geometry=geometry,
            )
        )

    if not features:
        raise ValueError(
            "No usable geometries found in the input vector file."
        )

    return features


def find_coord_name(
    ds: xr.Dataset,
    candidates: tuple[str, ...],
) -> str:
    """Find the first matching coordinate/variable name."""

    available = list(ds.coords) + list(ds.variables)

    lower_map = {
        name.lower(): name
        for name in available
    }

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    raise KeyError(
        f"Could not find a coordinate among {candidates}. "
        f"Available coordinates/variables: {available}"
    )


def find_time_name(
    ds: xr.Dataset,
) -> str | None:
    """Return the dataset time coordinate name if present."""

    for candidate in TIME_CANDIDATES:
        for name in list(ds.coords) + list(ds.variables):
            if name.lower() == candidate:
                return name

    return None


def choose_data_variable(
    ds: xr.Dataset,
    preferred: str,
) -> str:
    """
    Identify the requested climate variable while ignoring metadata
    variables such as CRS and time bounds.
    """

    for name in ds.data_vars:
        if name.lower() == preferred.lower():
            return name

    non_metadata_variables = [
        name
        for name in ds.data_vars
        if "bound" not in name.lower()
        and not name.lower().endswith("_bnds")
        and name.lower() != "crs"
    ]

    if len(non_metadata_variables) == 1:
        return non_metadata_variables[0]

    raise KeyError(
        f"Could not uniquely identify '{preferred}'. "
        f"Data variables: {list(ds.data_vars)}"
    )


def subset_time(
    ds: xr.Dataset,
    start_year: int,
    end_year: int,
) -> xr.Dataset:
    """Subset a dataset to an inclusive year range."""

    time_name = find_time_name(ds)

    if time_name is None:
        return ds

    start = f"{start_year}-01-01"
    end = f"{end_year}-12-31"

    try:
        return ds.sel(
            {
                time_name: slice(
                    start,
                    end,
                )
            }
        )

    except Exception:

        years = ds[time_name].dt.year

        return ds.where(
            (years >= start_year)
            & (years <= end_year),
            drop=True,
        )


def _copy_grid_mapping_and_bounds(
    *,
    source_ds: xr.Dataset,
    result_ds: xr.Dataset,
    variable_name: str,
) -> xr.Dataset:
    """
    Preserve source CRS and time-bounds metadata where NCI supplies them.

    Some NCI datasets genuinely omit CRS metadata. In those cases this
    function deliberately leaves the local subset unchanged; downstream
    raster processing has its own NARCliM CRS fallback.
    """

    grid_mapping_name = (
        source_ds[variable_name]
        .attrs
        .get("grid_mapping")
    )

    if (
        grid_mapping_name
        and grid_mapping_name in source_ds.variables
    ):
        result_ds[grid_mapping_name] = (
            source_ds[grid_mapping_name]
        )

        result_ds[
            variable_name
        ].attrs["grid_mapping"] = (
            grid_mapping_name
        )

    elif "crs" in source_ds.variables:

        result_ds["crs"] = source_ds["crs"]

        result_ds[
            variable_name
        ].attrs["grid_mapping"] = "crs"

    time_name = find_time_name(source_ds)

    if (
        time_name is not None
        and time_name in result_ds.coords
    ):
        bounds_name = (
            source_ds[time_name]
            .attrs
            .get("bounds")
        )

        if (
            bounds_name
            and bounds_name in source_ds.variables
        ):
            try:
                result_ds[bounds_name] = (
                    source_ds[bounds_name].sel(
                        {
                            time_name:
                                result_ds[time_name]
                        }
                    )
                )

            except Exception:
                result_ds[bounds_name] = (
                    source_ds[bounds_name]
                )

            result_ds[
                time_name
            ].attrs["bounds"] = bounds_name

    return result_ds


def _centres_to_corners(
    values: np.ndarray,
) -> np.ndarray:
    """
    Estimate grid-cell corner coordinates from a 2-D array of cell centres.

    NARCliM provides geographic longitude/latitude at cell centres.
    Cell polygons are needed so that even a very small intersection between
    a study polygon and a climate cell selects the complete climate cell.
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
            "At least a 2 x 2 grid is required to derive cell corners."
        )

    corners = np.empty(
        (
            n_rows + 1,
            n_cols + 1,
        ),
        dtype=np.float64,
    )

    # Interior corners: mean of the four neighbouring cell centres.
    corners[
        1:n_rows,
        1:n_cols,
    ] = (
        values[:-1, :-1]
        + values[:-1, 1:]
        + values[1:, :-1]
        + values[1:, 1:]
    ) / 4.0

    # Top and bottom edges.
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

    # Left and right edges.
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

    # Four outer corners.
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


def _candidate_window(
    *,
    lon: np.ndarray,
    lat: np.ndarray,
    geometry,
    padding_cells: int = 2,
) -> tuple[int, int, int, int]:
    """
    Find a conservative row/column window around a geographic geometry.

    This avoids testing every cell in the full NARCliM grid.
    """

    west, south, east, north = geometry.bounds

    # Estimate representative geographic centre-to-centre spacing.
    lon_dx = np.abs(
        np.diff(
            lon,
            axis=1,
        )
    )

    lat_dy = np.abs(
        np.diff(
            lat,
            axis=0,
        )
    )

    lon_spacing = float(
        np.nanmedian(lon_dx)
    )

    lat_spacing = float(
        np.nanmedian(lat_dy)
    )

    lon_padding = max(
        lon_spacing * padding_cells,
        0.05,
    )

    lat_padding = max(
        lat_spacing * padding_cells,
        0.05,
    )

    candidate_mask = (
        np.isfinite(lon)
        & np.isfinite(lat)
        & (lon >= west - lon_padding)
        & (lon <= east + lon_padding)
        & (lat >= south - lat_padding)
        & (lat <= north + lat_padding)
    )

    candidate_rows, candidate_cols = np.where(
        candidate_mask
    )

    if candidate_rows.size == 0:
        raise ValueError(
            "No NARCliM climate cells are geographically close to the "
            "supplied geometry. Check whether the selected climate domain "
            "covers the study area."
        )

    n_rows, n_cols = lon.shape

    row_start = max(
        int(candidate_rows.min()) - 1,
        0,
    )

    row_stop = min(
        int(candidate_rows.max()) + 2,
        n_rows,
    )

    col_start = max(
        int(candidate_cols.min()) - 1,
        0,
    )

    col_stop = min(
        int(candidate_cols.max()) + 2,
        n_cols,
    )

    return (
        row_start,
        row_stop,
        col_start,
        col_stop,
    )


def _polygon_intersection_mask(
    ds: xr.Dataset,
    geometry,
) -> np.ndarray:
    """
    Select every NARCliM climate cell intersected by the supplied polygon.

    NCI's own 2-D geographic ``lon``/``lat`` arrays are the spatial truth.
    A cell is selected if ANY part of its polygon intersects/touches the
    study polygon. The full climate-cell value is retained.
    """

    if "lon" not in ds.variables:
        raise ValueError(
            "Dataset does not contain the 2-D 'lon' coordinate."
        )

    if "lat" not in ds.variables:
        raise ValueError(
            "Dataset does not contain the 2-D 'lat' coordinate."
        )

    lon = np.asarray(
        ds["lon"].values,
        dtype=np.float64,
    )

    lat = np.asarray(
        ds["lat"].values,
        dtype=np.float64,
    )

    if lon.ndim != 2 or lat.ndim != 2:
        raise ValueError(
            "NARCliM polygon extraction requires 2-D lon/lat arrays."
        )

    if lon.shape != lat.shape:
        raise ValueError(
            "Longitude and latitude arrays have different shapes."
        )

    n_rows, n_cols = lon.shape

    # Build true geographic cell corners.
    lon_corners = _centres_to_corners(
        lon
    )

    lat_corners = _centres_to_corners(
        lat
    )

    (
        row_start,
        row_stop,
        col_start,
        col_stop,
    ) = _candidate_window(
        lon=lon,
        lat=lat,
        geometry=geometry,
    )

    selected = np.zeros(
        (
            n_rows,
            n_cols,
        ),
        dtype=bool,
    )

    prepared_geometry = prep(
        geometry
    )

    for row in range(
        row_start,
        row_stop,
    ):
        for col in range(
            col_start,
            col_stop,
        ):
            cell_lon = [
                lon_corners[row, col],
                lon_corners[row, col + 1],
                lon_corners[row + 1, col + 1],
                lon_corners[row + 1, col],
            ]

            cell_lat = [
                lat_corners[row, col],
                lat_corners[row, col + 1],
                lat_corners[row + 1, col + 1],
                lat_corners[row + 1, col],
            ]

            if not (
                np.all(
                    np.isfinite(cell_lon)
                )
                and np.all(
                    np.isfinite(cell_lat)
                )
            ):
                continue

            cell_polygon = Polygon(
                zip(
                    cell_lon,
                    cell_lat,
                )
            )

            if cell_polygon.is_empty:
                continue

            # ``buffer(0)`` repairs occasional minor polygon validity issues
            # caused by numerical corner reconstruction.
            if not cell_polygon.is_valid:
                cell_polygon = cell_polygon.buffer(0)

            if (
                not cell_polygon.is_empty
                and prepared_geometry.intersects(
                    cell_polygon
                )
            ):
                selected[
                    row,
                    col,
                ] = True

    return selected


def _crop_to_selected_cells(
    ds: xr.Dataset,
    *,
    mask_values: np.ndarray,
    variable_name: str,
) -> xr.Dataset:
    """
    Crop to the smallest continuous NARCliM window containing all selected
    cells while masking non-intersected cells inside that rectangle.
    """

    selected_rows, selected_cols = np.where(
        mask_values
    )

    if selected_rows.size == 0:
        raise ValueError(
            "The geometry does not intersect any NARCliM climate cells."
        )

    row_start = int(
        selected_rows.min()
    )

    row_stop = (
        int(
            selected_rows.max()
        )
        + 1
    )

    col_start = int(
        selected_cols.min()
    )

    col_stop = (
        int(
            selected_cols.max()
        )
        + 1
    )

    # Determine native horizontal dimension names from the 2-D lon array.
    lon = ds["lon"]

    if lon.ndim != 2:
        raise ValueError(
            "Expected the NARCliM longitude coordinate to be 2-D."
        )

    y_dim, x_dim = lon.dims

    cropped_ds = ds.isel(
        {
            y_dim: slice(
                row_start,
                row_stop,
            ),
            x_dim: slice(
                col_start,
                col_stop,
            ),
        }
    )

    cropped_mask = mask_values[
        row_start:row_stop,
        col_start:col_stop,
    ]

    mask_da = xr.DataArray(
        cropped_mask,
        dims=(
            y_dim,
            x_dim,
        ),
        coords={
            y_dim:
                cropped_ds.coords[y_dim],
            x_dim:
                cropped_ds.coords[x_dim],
        },
    )

    result = (
        cropped_ds[
            variable_name
        ]
        .where(
            mask_da
        )
        .to_dataset(
            name=variable_name
        )
    )

    # Preserve true geographic coordinates plus native grid coordinates.
    for coordinate_name in (
        "lon",
        "lat",
        "rlon",
        "rlat",
        "height",
        "time",
    ):
        if coordinate_name in cropped_ds.coords:

            result = result.assign_coords(
                {
                    coordinate_name:
                        cropped_ds[
                            coordinate_name
                        ]
                }
            )

        elif coordinate_name in cropped_ds.variables:

            result[
                coordinate_name
            ] = cropped_ds[
                coordinate_name
            ]

    result = _copy_grid_mapping_and_bounds(
        source_ds=cropped_ds,
        result_ds=result,
        variable_name=variable_name,
    )

    result.attrs = dict(
        cropped_ds.attrs
    )

    return result


def extract_polygon(
    ds: xr.Dataset,
    geometry,
    variable_name: str,
    *,
    verbose: bool = True,
) -> xr.Dataset:
    """
    Extract every NARCliM climate cell intersected by a polygon.

    This uses NCI's native 2-D geographic longitude/latitude coordinates,
    not transformed rlon/rlat coordinates.

    Even if only a tiny part of a climate cell touches the polygon, the
    complete climate cell is retained. This is suitable for both the
    ~4 km SEAus domain and the ~18 km AUS-18 domain.
    """

    if geometry is None:
        raise ValueError(
            "Polygon geometry is None."
        )

    if geometry.is_empty:
        raise ValueError(
            "Polygon geometry is empty."
        )

    mask_values = (
        _polygon_intersection_mask(
            ds,
            geometry,
        )
    )

    selected_cell_count = int(
        np.count_nonzero(
            mask_values
        )
    )

    if selected_cell_count == 0:
        raise ValueError(
            "The polygon does not intersect any NARCliM climate cells."
        )

    result = _crop_to_selected_cells(
        ds,
        mask_values=mask_values,
        variable_name=variable_name,
    )

    result.attrs[
        "spatial_selection_method"
    ] = "geographic_cell_intersection"

    result.attrs[
        "selected_grid_cells"
    ] = selected_cell_count

    if verbose:

        lon = np.asarray(
            result["lon"].values,
            dtype=float,
        )

        lat = np.asarray(
            result["lat"].values,
            dtype=float,
        )

        print(
            "[SPATIAL METHOD] geographic_cell_intersection | "
            f"selected={selected_cell_count} climate cell(s)",
            flush=True,
        )

        print(
            "[SPATIAL EXTENT] "
            f"west={np.nanmin(lon):.6f}; "
            f"south={np.nanmin(lat):.6f}; "
            f"east={np.nanmax(lon):.6f}; "
            f"north={np.nanmax(lat):.6f}",
            flush=True,
        )

    return result


def extract_point(
    ds: xr.Dataset,
    geometry,
    variable_name: str,
    *,
    verbose: bool = True,
) -> xr.Dataset:
    """
    Extract the nearest NARCliM climate cell for a point.

    NCI's 2-D geographic lon/lat cell centres are used directly, avoiding
    any ambiguity in the native rotated-grid coordinates.
    """

    if geometry is None or geometry.is_empty:
        raise ValueError(
            "Point geometry is empty."
        )

    # Accept a true Point. For safety, use a representative point if a
    # point-like geometry reaches this function unexpectedly.
    if isinstance(
        geometry,
        Point,
    ):
        point = geometry
    else:
        point = geometry.representative_point()

    if (
        "lon" not in ds.variables
        or "lat" not in ds.variables
    ):
        raise ValueError(
            "Dataset must contain 2-D lon/lat coordinates for point extraction."
        )

    lon = np.asarray(
        ds["lon"].values,
        dtype=np.float64,
    )

    lat = np.asarray(
        ds["lat"].values,
        dtype=np.float64,
    )

    if lon.ndim != 2 or lat.ndim != 2:
        raise ValueError(
            "NARCliM point extraction requires 2-D lon/lat arrays."
        )

    target_lon = float(
        point.x
    )

    target_lat = float(
        point.y
    )

    # Approximate geographic distance. Longitude is scaled by latitude,
    # which is sufficient for choosing the nearest grid centre locally.
    longitude_scale = np.cos(
        np.deg2rad(
            target_lat
        )
    )

    distance_squared = (
        (
            (lon - target_lon)
            * longitude_scale
        ) ** 2
        + (lat - target_lat) ** 2
    )

    distance_squared[
        ~np.isfinite(
            distance_squared
        )
    ] = np.inf

    flat_index = int(
        np.argmin(
            distance_squared
        )
    )

    row, col = np.unravel_index(
        flat_index,
        distance_squared.shape,
    )

    y_dim, x_dim = ds["lon"].dims

    # Keep singleton spatial dimensions so downstream NetCDF structure
    # remains predictable.
    result = ds.isel(
        {
            y_dim: slice(
                row,
                row + 1,
            ),
            x_dim: slice(
                col,
                col + 1,
            ),
        }
    )

    # Keep only the requested climate variable plus required coordinates
    # and metadata.
    output = result[
        variable_name
    ].to_dataset(
        name=variable_name
    )

    for coordinate_name in (
        "lon",
        "lat",
        "rlon",
        "rlat",
        "height",
        "time",
    ):
        if coordinate_name in result.coords:

            output = output.assign_coords(
                {
                    coordinate_name:
                        result[
                            coordinate_name
                        ]
                }
            )

        elif coordinate_name in result.variables:

            output[
                coordinate_name
            ] = result[
                coordinate_name
            ]

    output = _copy_grid_mapping_and_bounds(
        source_ds=result,
        result_ds=output,
        variable_name=variable_name,
    )

    output.attrs = dict(
        result.attrs
    )

    output.attrs[
        "spatial_selection_method"
    ] = "nearest_geographic_cell"

    output.attrs[
        "selected_grid_cells"
    ] = 1

    if verbose:

        selected_lon = float(
            np.asarray(
                result["lon"].values
            ).squeeze()
        )

        selected_lat = float(
            np.asarray(
                result["lat"].values
            ).squeeze()
        )

        print(
            "[SPATIAL METHOD] nearest_geographic_cell | "
            f"selected_lon={selected_lon:.6f}; "
            f"selected_lat={selected_lat:.6f}",
            flush=True,
        )

    return output