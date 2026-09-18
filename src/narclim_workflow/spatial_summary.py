"""
spatial_summary.py

General spatial summarisation of NARCliM ensemble hazard rasters for user-supplied
POLYGON or POINT feature datasets.

The module reads the original pixel-wise ensemble rasters created by
``uncertainty.py`` from:

    Climate_Indicies/
        Ensemble_Stats/

For every available combination of:

    variable × scenario × time_horizon × source feature

the uncertainty workflow provides three rasters:

    <variable>_Min_<scenario>_<start>_<end>.tif
    <variable>_Mean_<scenario>_<start>_<end>.tif
    <variable>_Max_<scenario>_<start>_<end>.tif

These rasters already represent climate-model uncertainty AFTER temporal
averaging within each requested climate horizon.

POLYGON INPUT
-------------
Every raster pixel touched by the polygon is included (``all_touched=True``).

For each polygon:

    min  = spatial mean of valid cells in the ensemble-min raster
    mean = spatial mean of valid cells in the ensemble-mean raster
    max  = spatial mean of valid cells in the ensemble-max raster

POINT INPUT
-----------
For each point, the module samples the raster pixel containing that point:

    min  = value from ensemble-min raster
    mean = value from ensemble-mean raster
    max  = value from ensemble-max raster

OUTPUT
------
A GeoPackage is written to:

    Climate_Indicies/
        Hazards/
            <input_dataset_name>_Hazards.gpkg

The GeoPackage contains ONE LAYER PER CLIMATE VARIABLE.

Each output record also contains ``time_horizon_year``, the representative
year of the climate period. For example, 1985-2014 is represented by 2000,
2041-2060 by 2050, and 2080-2100 by 2090.

A non-spatial CSV containing all variable records is also written beside the
GeoPackage using the same base name, for example:

    Planning_Zones_Hazards.csv
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask


SPATIAL_SUMMARY_CODE_VERSION = "2026-09-18-optional-mean-change-fix-v7"

SUPPORTED_GEOMETRIES = {
    "Point",
    "MultiPoint",
    "Polygon",
    "MultiPolygon",
}

STATISTICS = ("min", "mean", "max")


@dataclass(frozen=True)
class EnsembleRasterGroup:
    """Three ensemble-statistic rasters belonging to one climate combination."""

    variable: str
    scenario: str
    time_horizon: str
    source_feature_id: str
    min_path: Path
    mean_path: Path | None
    max_path: Path


def _log(message: str, *, verbose: bool) -> None:
    if verbose:
        print(message, flush=True)


def _safe_name(value: object, *, max_length: int = 60) -> str:
    text = str(value).strip()
    text = re.sub(r"[^0-9A-Za-z_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")

    if not text:
        text = "layer"

    if text[0].isdigit():
        text = f"v_{text}"

    return text[:max_length]



_STORM_DERIVED_RE = re.compile(
    r"^StormDaysGT\d+(?:p\d+)?$",
    flags=re.IGNORECASE,
)


def _output_decimals(
    variable: str,
) -> int:
    """
    Use 3 decimals for rare storm-day frequencies and 1 decimal for
    other climate variables.
    """

    if _STORM_DERIVED_RE.match(
        str(variable)
    ):
        return 3

    return 1



def _representative_year(
    start_year: int,
    end_year: int,
) -> int:
    """
    Return the representative midpoint year for a climate time horizon.

    Examples
    --------
    1985-2014 -> 2000
    2041-2060 -> 2050
    2040-2060 -> 2050
    2080-2100 -> 2090
    """

    return int(
        round(
            (
                int(start_year)
                + int(end_year)
            )
            / 2.0
        )
    )


def _extract_period_from_raster_name(
    raster_path: Path,
) -> tuple[int, int] | None:
    """
    Extract start/end years from the new ensemble raster filename.

    Example
    -------
    TXge35_Mean_ssp245_2040_2060.tif -> (2040, 2060)
    """

    match = re.search(
        r"_(?:historical|ssp[0-9]+)_"
        r"(?P<start_year>[0-9]{4})_"
        r"(?P<end_year>[0-9]{4})"
        r"(?:_|$)",
        raster_path.stem,
        flags=re.IGNORECASE,
    )

    if match is None:
        return None

    return (
        int(match.group("start_year")),
        int(match.group("end_year")),
    )


def _get_time_horizon_year(
    *,
    time_horizon: str,
    time_windows: dict[str, tuple[int, int]] | None,
    raster_path: Path,
) -> int | None:
    """
    Return the representative year for one climate time horizon.

    Priority
    --------
    1. Use ``time_windows`` when supplied.
    2. Otherwise derive the period from the ensemble raster filename.
    """

    if (
        time_windows is not None
        and time_horizon in time_windows
    ):
        start_year, end_year = time_windows[
            time_horizon
        ]

        return _representative_year(
            start_year,
            end_year,
        )

    period = _extract_period_from_raster_name(
        raster_path
    )

    if period is None:
        return None

    return _representative_year(
        period[0],
        period[1],
    )


def _format_time_horizon(
    time_horizon: str,
    time_windows: dict[str, tuple[int, int]] | None,
    raster_path: Path | None = None,
) -> str:
    """
    Format the time-horizon label and append its representative year.

    Examples
    --------
    baseline -> baseline (2000)
    mid_term -> mid_term (2050)
    long_term -> long_term (2090)
    """

    representative_year = None

    if (
        time_windows is not None
        and time_horizon in time_windows
    ):
        start_year, end_year = time_windows[
            time_horizon
        ]

        representative_year = _representative_year(
            start_year,
            end_year,
        )

    elif raster_path is not None:

        period = _extract_period_from_raster_name(
            raster_path
        )

        if period is not None:
            representative_year = _representative_year(
                period[0],
                period[1],
            )

    if representative_year is None:
        return time_horizon

    return (
        f"{time_horizon} "
        f"({representative_year})"
    )

def _read_input_features(
    *,
    feature_path: Path,
    feature_layer: str | None,
    feature_id_field: str | None,
) -> tuple[gpd.GeoDataFrame, str, str]:
    if not feature_path.exists():
        raise FileNotFoundError(
            f"Input feature dataset was not found:\n{feature_path}"
        )

    if feature_layer is None:
        features = gpd.read_file(feature_path)
    else:
        features = gpd.read_file(feature_path, layer=feature_layer)

    if features.empty:
        raise ValueError("The input feature dataset contains no features.")

    if features.crs is None:
        raise ValueError("The input feature dataset has no CRS.")

    geometry_types = set(
        features.geometry.geom_type.dropna().unique()
    )

    unsupported = geometry_types - SUPPORTED_GEOMETRIES

    if unsupported:
        raise ValueError(
            "Only point and polygon feature inputs are supported. "
            f"Unsupported geometry type(s): {sorted(unsupported)}"
        )

    has_points = bool(geometry_types & {"Point", "MultiPoint"})
    has_polygons = bool(geometry_types & {"Polygon", "MultiPolygon"})

    if has_points and has_polygons:
        raise ValueError(
            "The input layer mixes point and polygon geometries. "
            "Please provide a layer containing only one geometry family."
        )

    if has_polygons:
        geometry_mode = "polygon"

        invalid = ~features.geometry.is_valid

        if invalid.any():
            try:
                features.loc[invalid, "geometry"] = (
                    features.loc[invalid, "geometry"].make_valid()
                )
            except AttributeError:
                features.loc[invalid, "geometry"] = (
                    features.loc[invalid, "geometry"].buffer(0)
                )

    elif has_points:
        geometry_mode = "point"

    else:
        raise ValueError(
            "No usable point or polygon geometries were found."
        )

    features = features.copy()

    if feature_id_field is not None:
        if feature_id_field not in features.columns:
            raise ValueError(
                f"feature_id_field={feature_id_field!r} was not found.\n"
                f"Available fields: {list(features.columns)}"
            )

        if features[feature_id_field].isna().any():
            raise ValueError(
                f"feature_id_field={feature_id_field!r} contains null values."
            )
    else:
        feature_id_field = "feature_id"

        if feature_id_field in features.columns:
            i = 1
            while f"feature_id_{i}" in features.columns:
                i += 1
            feature_id_field = f"feature_id_{i}"

        features[feature_id_field] = np.arange(
            1,
            len(features) + 1,
            dtype=np.int64,
        )

    return features, feature_id_field, geometry_mode


def _parse_ensemble_raster_name(
    raster_path: Path,
) -> tuple[str, str, str, str, str] | None:
    """
    Parse an ensemble-statistic raster filename.

    Supported naming conventions
    ----------------------------
    NEW:
        TXge35_Min_historical_1985_2014.tif
        TXge35_Mean_ssp245_2040_2060.tif
        TXge35_Max_ssp370_2080_2100.tif

    OLD:
        TXge35_ssp245_mid_term_ID1_ensemble_mean.tif

    Returns
    -------
    tuple or None
        variable, scenario, time_horizon, source_feature_id, statistic

    Notes
    -----
    For the new filename convention, the time-horizon name is read from
    the GeoTIFF metadata written by uncertainty.py. If it is unavailable,
    a year-range label is used as a safe fallback.

    The normal single-feature workflow omits ID1 from the new filename, so
    source_feature_id defaults to "ID1".
    """

    stem = raster_path.stem

    # ==================================================================
    # 1. NEW NAMING CONVENTION
    # ==================================================================

    new_match = re.match(
        r"^(?P<variable>.+?)_"
        r"(?P<stat>Min|Mean|Max)_"
        r"(?P<scenario>historical|ssp[0-9]+)_"
        r"(?P<start_year>[0-9]{4})_"
        r"(?P<end_year>[0-9]{4})"
        r"(?:_(?P<feature_id>.+))?$",
        stem,
        flags=re.IGNORECASE,
    )

    if new_match is not None:

        variable = new_match.group("variable")
        scenario = new_match.group("scenario")
        statistic = new_match.group("stat").lower()
        start_year = new_match.group("start_year")
        end_year = new_match.group("end_year")

        source_feature_id = (
            new_match.group("feature_id")
            or "ID1"
        )

        # Read the original horizon name from the raster tags.
        # uncertainty.py writes this metadata when creating the raster.
        time_horizon = None

        try:
            with rasterio.open(raster_path) as src:
                tags = src.tags()
                time_horizon = tags.get(
                    "time_horizon"
                )
        except Exception:
            time_horizon = None

        if not time_horizon:
            time_horizon = (
                f"{start_year}_{end_year}"
            )

        return (
            variable,
            scenario,
            time_horizon,
            source_feature_id,
            statistic,
        )

    # ==================================================================
    # 2. OLD NAMING CONVENTION
    # ==================================================================

    statistic_match = re.match(
        r"^(?P<body>.+)_ensemble_(?P<stat>min|mean|max)$",
        stem,
        flags=re.IGNORECASE,
    )

    if statistic_match is None:
        return None

    body = statistic_match.group("body")
    statistic = (
        statistic_match.group("stat")
        .lower()
    )

    scenario_match = re.match(
        r"^(?P<variable>.+?)_"
        r"(?P<scenario>historical|ssp[0-9]+)_"
        r"(?P<rest>.+)$",
        body,
        flags=re.IGNORECASE,
    )

    if scenario_match is None:
        return None

    variable = scenario_match.group(
        "variable"
    )

    scenario = scenario_match.group(
        "scenario"
    )

    rest = scenario_match.group(
        "rest"
    )

    if "_" not in rest:
        return None

    time_horizon, source_feature_id = (
        rest.rsplit(
            "_",
            1,
        )
    )

    return (
        variable,
        scenario,
        time_horizon,
        source_feature_id,
        statistic,
    )

def _discover_ensemble_rasters(
    *,
    ensemble_folder: Path,
    variables: list[str] | None,
    scenarios: list[str] | None,
    time_horizons: list[str] | None,
    source_feature_ids: list[str] | None,
    include_mean: bool = True,
) -> list[EnsembleRasterGroup]:
    if not ensemble_folder.exists():
        raise FileNotFoundError(
            f"Ensemble raster folder was not found:\n{ensemble_folder}"
        )

    grouped: dict[
        tuple[str, str, str, str],
        dict[str, Path],
    ] = {}

    for raster_path in sorted(ensemble_folder.glob("*.tif")):
        parsed = _parse_ensemble_raster_name(raster_path)

        if parsed is None:
            continue

        variable, scenario, time_horizon, source_feature_id, statistic = parsed

        # uncertainty.py may also write pre-calculated *_Change rasters.
        # spatial_summary.py calculates change from the original climate
        # variable rasters itself, so do not treat those change rasters as
        # separate climate variables here.
        if variable.endswith("_Change"):
            continue

        if variables is not None and variable not in variables:
            continue

        if scenarios is not None and scenario not in scenarios:
            continue

        if time_horizons is not None and time_horizon not in time_horizons:
            continue

        if (
            source_feature_ids is not None
            and source_feature_id not in source_feature_ids
        ):
            continue

        key = (
            variable,
            scenario,
            time_horizon,
            source_feature_id,
        )

        grouped.setdefault(key, {})[statistic] = raster_path

    complete_groups: list[EnsembleRasterGroup] = []

    for (
        variable,
        scenario,
        time_horizon,
        source_feature_id,
    ), paths in sorted(grouped.items()):

        required_statistics = {"min", "max"}
        if include_mean:
            required_statistics.add("mean")

        missing = required_statistics - set(paths)

        if missing:
            continue

        complete_groups.append(
            EnsembleRasterGroup(
                variable=variable,
                scenario=scenario,
                time_horizon=time_horizon,
                source_feature_id=source_feature_id,
                min_path=paths["min"],
                mean_path=paths.get("mean"),
                max_path=paths["max"],
            )
        )

    if not complete_groups:
        raise ValueError(
            "No complete ensemble raster groups were found for the requested "
            "filters. Min and Max are required; Mean is required only when "
            "include_mean=True."
        )

    return complete_groups


def _validate_raster_group(
    group: EnsembleRasterGroup,
) -> tuple[object, object, tuple[int, int]]:
    reference_crs = None
    reference_transform = None
    reference_shape = None

    raster_paths = [group.min_path, group.max_path]
    if group.mean_path is not None:
        raster_paths.insert(1, group.mean_path)

    for path in raster_paths:
        with rasterio.open(path) as src:
            if src.crs is None:
                raise ValueError(f"Raster has no CRS:\n{path}")

            shape = (src.height, src.width)

            if reference_crs is None:
                reference_crs = src.crs
                reference_transform = src.transform
                reference_shape = shape
            else:
                if src.crs != reference_crs:
                    raise ValueError(
                        "The min/mean/max ensemble rasters have different CRS "
                        f"for {group.variable} | {group.scenario} | "
                        f"{group.time_horizon}."
                    )

                if not src.transform.almost_equals(reference_transform):
                    raise ValueError(
                        "The min/mean/max ensemble rasters have different grid "
                        f"alignment for {group.variable} | {group.scenario} | "
                        f"{group.time_horizon}."
                    )

                if shape != reference_shape:
                    raise ValueError(
                        "The min/mean/max ensemble rasters have different shapes "
                        f"for {group.variable} | {group.scenario} | "
                        f"{group.time_horizon}."
                    )

    return reference_crs, reference_transform, reference_shape


def _polygon_spatial_mean(
    *,
    raster_path: Path,
    geometry,
    all_touched: bool,
) -> tuple[float, int]:
    with rasterio.open(raster_path) as src:
        raster = src.read(1, masked=True)

        if geometry is None or geometry.is_empty:
            return np.nan, 0

        inside = geometry_mask(
            [geometry],
            out_shape=(src.height, src.width),
            transform=src.transform,
            invert=True,
            all_touched=all_touched,
        )

        combined_mask = (
            np.ma.getmaskarray(raster)
            | ~inside
        )

        values = np.ma.array(
            raster,
            mask=combined_mask,
        ).compressed()

        if values.size == 0:
            return np.nan, 0

        return float(np.nanmean(values)), int(values.size)


def _sample_point_value(
    *,
    raster_path: Path,
    point,
) -> float:
    if point is None or point.is_empty:
        return np.nan

    if point.geom_type == "MultiPoint":
        point = list(point.geoms)[0]

    x = float(point.x)
    y = float(point.y)

    with rasterio.open(raster_path) as src:
        bounds = src.bounds

        if not (
            bounds.left <= x <= bounds.right
            and bounds.bottom <= y <= bounds.top
        ):
            return np.nan

        sampled = next(
            src.sample(
                [(x, y)],
                masked=True,
            )
        )

        value = sampled[0]

        if np.ma.is_masked(value):
            return np.nan

        value = float(value)

        if (
            src.nodata is not None
            and np.isclose(value, src.nodata)
        ):
            return np.nan

        if not np.isfinite(value):
            return np.nan

        return value


def _summarise_group_for_features(
    *,
    features_original: gpd.GeoDataFrame,
    geometry_mode: str,
    group: EnsembleRasterGroup,
    keep_fields: list[str],
    time_windows: dict[str, tuple[int, int]] | None,
    all_touched: bool,
    include_mean: bool = True,
) -> gpd.GeoDataFrame:
    raster_crs, _, _ = _validate_raster_group(group)

    features_raster = features_original.to_crs(raster_crs)

    output_records = []

    formatted_horizon = _format_time_horizon(
        group.time_horizon,
        time_windows,
        raster_path=group.min_path,
    )

    # Representative year for direct GIS/chart use.
    #
    # Examples:
    #   1985-2014 -> 2000
    #   2041-2060 -> 2050
    #   2080-2100 -> 2090
    time_horizon_year = _get_time_horizon_year(
        time_horizon=group.time_horizon,
        time_windows=time_windows,
        raster_path=group.min_path,
    )

    for position in range(len(features_original)):
        original_row = features_original.iloc[position]
        raster_row = features_raster.iloc[position]
        geometry = raster_row.geometry

        if geometry_mode == "polygon":
            min_value, min_count = _polygon_spatial_mean(
                raster_path=group.min_path,
                geometry=geometry,
                all_touched=all_touched,
            )

            if include_mean:
                mean_value, mean_count = _polygon_spatial_mean(
                    raster_path=group.mean_path,
                    geometry=geometry,
                    all_touched=all_touched,
                )
            else:
                mean_value, mean_count = np.nan, 0

            max_value, max_count = _polygon_spatial_mean(
                raster_path=group.max_path,
                geometry=geometry,
                all_touched=all_touched,
            )

            n_pixels = int(
                max(
                    min_count,
                    mean_count,
                    max_count,
                )
            )
        else:
            min_value = _sample_point_value(
                raster_path=group.min_path,
                point=geometry,
            )

            if include_mean:
                mean_value = _sample_point_value(
                    raster_path=group.mean_path,
                    point=geometry,
                )
            else:
                mean_value = np.nan

            max_value = _sample_point_value(
                raster_path=group.max_path,
                point=geometry,
            )

            n_pixels = int(
                any(
                    np.isfinite(value)
                    for value in (
                        min_value,
                        mean_value,
                        max_value,
                    )
                )
            )

        # --------------------------------------------------------------
        # Final output precision:
        #   most variables -> 1 decimal
        #   StormDaysGT*    -> 3 decimals
        #
        # This prevents rare storm frequencies such as 0.033 days/year
        # from being rounded to zero.
        # --------------------------------------------------------------
        output_decimals = _output_decimals(
            group.variable
        )

        min_value = (
            round(
                float(min_value),
                output_decimals,
            )
            if np.isfinite(min_value)
            else np.nan
        )

        mean_value = (
            round(
                float(mean_value),
                output_decimals,
            )
            if np.isfinite(mean_value)
            else np.nan
        )

        max_value = (
            round(
                float(max_value),
                output_decimals,
            )
            if np.isfinite(max_value)
            else np.nan
        )

        record = {
            field: original_row[field]
            for field in keep_fields
        }

        record.update(
            {
                "variable": group.variable,
                "scenario": group.scenario,
                "time_horizon": formatted_horizon,
                "time_horizon_year": time_horizon_year,
                "source_feature": group.source_feature_id,
                "min": min_value,
                **({"mean": mean_value} if include_mean else {}),
                "max": max_value,
                "n_pixels": n_pixels,
                "geometry": original_row.geometry,
            }
        )

        if include_mean:
            record["mean"] = mean_value

        output_records.append(record)

    return gpd.GeoDataFrame(
        output_records,
        geometry="geometry",
        crs=features_original.crs,
    )


def run_spatial_summary(
    *,
    feature_path: str | Path,
    climate_indices_root: str | Path,
    feature_id_field: str | None = None,
    feature_layer: str | None = None,
    dataset_name: str | None = None,
    output_name: str | None = None,
    variables: list[str] | None = None,
    scenarios: list[str] | None = None,
    time_horizons: list[str] | None = None,
    source_feature_ids: list[str] | None = None,
    time_windows: dict[str, tuple[int, int]] | None = None,
    keep_fields: list[str] | None = None,
    all_touched: bool = True,
    include_mean: bool = True,
    calculate_change: bool = True,
    baseline_scenario: str = "historical",
    baseline_time_horizon: str = "baseline",
    write_feature_summary_layer: bool = True,
    overwrite: bool = True,
    verbose: bool = True,
) -> dict[str, Path]:
    """
    Summarise NARCliM ensemble hazard rasters for polygon OR point features.

    Polygon inputs:
        - use original ensemble rasters directly
        - include all intersected cells when all_touched=True
        - calculate spatial mean for ensemble min / mean / max rasters

    Point inputs:
        - sample the raster cell containing each point
        - assign min / mean / max directly

    Optional additions:
        - include_mean=False omits ensemble-mean values from the outputs.
        - calculate_change=True calculates scenario minus historical-baseline
          change for min/max and, when requested, mean.
        - write_feature_summary_layer=True writes a wide copy of FEATURE_PATH
          with the calculated climate-change fields. The input file is never modified.
    """

    feature_path = Path(feature_path)
    climate_indices_root = Path(climate_indices_root)

    ensemble_folder = (
        climate_indices_root
        / "Ensemble_Stats"
    )

    hazards_folder = (
        climate_indices_root
        / "Hazards"
    )

    hazards_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    _log("=" * 100, verbose=verbose)
    _log("NARCliM SPATIAL SUMMARY WORKFLOW", verbose=verbose)
    _log("=" * 100, verbose=verbose)
    _log(
        f"[CODE VERSION] {SPATIAL_SUMMARY_CODE_VERSION}",
        verbose=verbose,
    )

    features, feature_id_field, geometry_mode = _read_input_features(
        feature_path=feature_path,
        feature_layer=feature_layer,
        feature_id_field=feature_id_field,
    )

    _log(f"[INPUT] {feature_path}", verbose=verbose)
    _log(f"[FEATURES] {len(features)}", verbose=verbose)
    _log(
        f"[GEOMETRY MODE] {geometry_mode.upper()}",
        verbose=verbose,
    )
    _log(f"[INPUT CRS] {features.crs}", verbose=verbose)

    if geometry_mode == "polygon":
        _log(
            f"[POLYGON METHOD] raster zonal mean | "
            f"all_touched={all_touched}",
            verbose=verbose,
        )
    else:
        _log(
            "[POINT METHOD] sample raster pixel containing each point",
            verbose=verbose,
        )

    non_geometry_fields = [
        field
        for field in features.columns
        if field != "geometry"
    ]

    if keep_fields is None:
        keep_fields_final = non_geometry_fields
    else:
        missing = [
            field
            for field in keep_fields
            if field not in features.columns
        ]

        if missing:
            raise ValueError(
                "Unknown keep_fields: "
                + ", ".join(missing)
            )

        keep_fields_final = list(
            dict.fromkeys(
                [
                    feature_id_field,
                    *keep_fields,
                ]
            )
        )

    groups = _discover_ensemble_rasters(
        ensemble_folder=ensemble_folder,
        variables=variables,
        scenarios=scenarios,
        time_horizons=time_horizons,
        source_feature_ids=source_feature_ids,
        include_mean=include_mean,
    )

    _log(
        f"[ENSEMBLE GROUPS] {len(groups)}",
        verbose=verbose,
    )

    _log(
        "[VARIABLES] "
        + ", ".join(
            sorted(
                {
                    group.variable
                    for group in groups
                }
            )
        ),
        verbose=verbose,
    )

    if dataset_name is None:
        dataset_name = feature_path.stem

    if output_name is None:
        output_name = (
            f"{_safe_name(dataset_name)}"
            "_Hazards.gpkg"
        )

    if not output_name.lower().endswith(".gpkg"):
        output_name += ".gpkg"

    output_gpkg = (
        hazards_folder
        / output_name
    )

    if output_gpkg.exists():
        if overwrite:
            output_gpkg.unlink()
        else:
            raise FileExistsError(
                f"Output exists:\n{output_gpkg}\n"
                "Set overwrite=True to replace it."
            )

    results_by_variable: dict[
        str,
        list[gpd.GeoDataFrame],
    ] = {}

    for group_index, group in enumerate(
        groups,
        start=1,
    ):
        _log("-" * 100, verbose=verbose)

        _log(
            f"[GROUP {group_index}/{len(groups)}] "
            f"{group.variable} | "
            f"{group.scenario} | "
            f"{group.time_horizon} | "
            f"{group.source_feature_id}",
            verbose=verbose,
        )

        summary = _summarise_group_for_features(
            features_original=features,
            geometry_mode=geometry_mode,
            group=group,
            keep_fields=keep_fields_final,
            time_windows=time_windows,
            all_touched=all_touched,
            include_mean=include_mean,
        )

        results_by_variable.setdefault(
            group.variable,
            [],
        ).append(summary)

        valid_field = "mean" if include_mean else "min"
        valid_count = int(
            summary[valid_field].notna().sum()
        )

        _log(
            f"[OK] {valid_count}/{len(summary)} feature record(s) "
            "have climate values.",
            verbose=verbose,
        )

    # ==================================================================
    # BASELINE SUMMARIES FOR CHANGE CALCULATION
    # ==================================================================
    # Change is defined as:
    #     future scenario value - historical baseline value
    #
    # Baseline groups are discovered separately so users can request only
    # ssp245/ssp370 and mid/long-term outputs while change calculations still
    # have access to the historical baseline.
    baseline_by_variable: dict[str, list[gpd.GeoDataFrame]] = {}

    if calculate_change:
        baseline_groups = _discover_ensemble_rasters(
            ensemble_folder=ensemble_folder,
            variables=variables,
            scenarios=[baseline_scenario],
            time_horizons=[baseline_time_horizon],
            source_feature_ids=source_feature_ids,
            include_mean=include_mean,
        )

        for baseline_group in baseline_groups:
            baseline_summary = _summarise_group_for_features(
                features_original=features,
                geometry_mode=geometry_mode,
                group=baseline_group,
                keep_fields=keep_fields_final,
                time_windows=time_windows,
                all_touched=all_touched,
                include_mean=include_mean,
            )
            baseline_by_variable.setdefault(
                baseline_group.variable, []
            ).append(baseline_summary)

    written: dict[str, Path] = {}

    # Collect non-spatial copies of every variable layer so a single CSV
    # can be written alongside the multi-layer GeoPackage.
    csv_frames: list[pd.DataFrame] = []
    feature_change_frames: list[pd.DataFrame] = []

    for variable in sorted(results_by_variable):
        frames = results_by_variable[variable]

        combined = pd.concat(
            frames,
            ignore_index=True,
        )

        combined = gpd.GeoDataFrame(
            combined,
            geometry="geometry",
            crs=features.crs,
        )

        # --------------------------------------------------------------
        # Calculate change relative to the historical baseline.
        # Matching is by input feature ID and source feature ID.
        # --------------------------------------------------------------
        if calculate_change:
            baseline_frames = baseline_by_variable.get(variable, [])

            if not baseline_frames:
                _log(
                    f"[CHANGE SKIP] No baseline summary for {variable!r}. "
                    f"Expected {baseline_scenario!r} / "
                    f"{baseline_time_horizon!r}; the normal spatial summary "
                    "will still be written, but change fields cannot be "
                    "calculated for this variable.",
                    verbose=verbose,
                )
            else:
                baseline_table = pd.concat(
                    baseline_frames, ignore_index=True
                )

                key_fields = [feature_id_field, "source_feature"]
                stat_fields = ["min", "max"]
                if include_mean:
                    stat_fields.insert(1, "mean")

                baseline_values = baseline_table[
                    [*key_fields, *stat_fields]
                ].copy()
                baseline_values = baseline_values.drop_duplicates(
                    subset=key_fields
                )
                baseline_values = baseline_values.rename(
                    columns={
                        field: f"baseline_{field}"
                        for field in stat_fields
                    }
                )

                combined = combined.merge(
                    baseline_values,
                    on=key_fields,
                    how="left",
                )
                combined = gpd.GeoDataFrame(
                    combined, geometry="geometry", crs=features.crs
                )

                for field in stat_fields:
                    combined[f"change_{field}"] = (
                        pd.to_numeric(combined[field], errors="coerce")
                        - pd.to_numeric(
                            combined[f"baseline_{field}"], errors="coerce"
                        )
                    )


        # --------------------------------------------------------------
        # Keep GeoPackage/CSV precision consistent with the variable.
        # --------------------------------------------------------------
        output_decimals = _output_decimals(
            variable
        )

        climate_fields = ["min", "max"]
        if include_mean:
            climate_fields.insert(1, "mean")

        round_fields = list(climate_fields)
        if calculate_change:
            round_fields.extend(
                [f"baseline_{field}" for field in climate_fields]
            )
            round_fields.extend(
                [f"change_{field}" for field in climate_fields]
            )

        for climate_field in round_fields:
            if climate_field in combined.columns:
                combined[
                    climate_field
                ] = pd.to_numeric(
                    combined[
                        climate_field
                    ],
                    errors="coerce",
                ).round(
                    output_decimals
                )

        # Build a non-spatial copy for the combined CSV.
        csv_frame = pd.DataFrame(
            combined.drop(
                columns="geometry"
            )
        )

        csv_frames.append(
            csv_frame
        )

        if calculate_change and write_feature_summary_layer:
            change_fields = [
                field for field in combined.columns
                if field.startswith("change_")
            ]
            for _, change_row in combined.iterrows():
                # Baseline records are useful in the variable layers but do not
                # need duplicate zero-change columns in the wide feature layer.
                if str(change_row["scenario"]).lower() == str(baseline_scenario).lower():
                    continue

                row_data = {feature_id_field: change_row[feature_id_field]}
                for field in change_fields:
                    field_name = _safe_name(
                        f"{variable}_{change_row['scenario']}_"
                        f"{change_row['time_horizon_year']}_"
                        f"{change_row['source_feature']}_{field}",
                        max_length=60,
                    )
                    row_data[field_name] = change_row[field]
                feature_change_frames.append(pd.DataFrame([row_data]))

        layer_name = _safe_name(variable)

        combined.to_file(
            output_gpkg,
            layer=layer_name,
            driver="GPKG",
            index=False,
        )

        written[variable] = output_gpkg

        _log(
            f"[WRITE] layer={layer_name} | "
            f"records={len(combined)}",
            verbose=verbose,
        )

    # ==================================================================
    # WRITE FEATURE-PATH CHANGE SUMMARY LAYER
    # ==================================================================
    # This is a copy of the supplied FEATURE_PATH geometry/attributes with
    # climate-change fields joined to it. The original input dataset is not
    # edited.
    if calculate_change and write_feature_summary_layer and feature_change_frames:
        feature_changes_long = pd.concat(
            feature_change_frames, ignore_index=True
        )
        feature_changes_wide = feature_changes_long.groupby(
            feature_id_field, as_index=False
        ).first()

        feature_summary = features.merge(
            feature_changes_wide,
            on=feature_id_field,
            how="left",
        )
        feature_summary = gpd.GeoDataFrame(
            feature_summary, geometry="geometry", crs=features.crs
        )

        feature_layer_name = _safe_name(
            f"{dataset_name}_Change_Summary"
        )
        feature_summary.to_file(
            output_gpkg,
            layer=feature_layer_name,
            driver="GPKG",
            index=False,
        )
        written["feature_change_summary"] = output_gpkg

        _log(
            f"[WRITE] layer={feature_layer_name} | "
            f"records={len(feature_summary)}",
            verbose=verbose,
        )

    # ==================================================================
    # WRITE COMBINED CSV
    #
    # Example:
    #   Planning_Zones_Hazards.gpkg
    #   Planning_Zones_Hazards.csv
    #
    # The GeoPackage keeps one spatial layer per climate variable.
    # The CSV combines all climate-variable records into one table and
    # excludes geometry because CSV is non-spatial.
    # ==================================================================

    output_csv = (
        output_gpkg.with_suffix(
            ".csv"
        )
    )

    if csv_frames:

        combined_csv = pd.concat(
            csv_frames,
            ignore_index=True,
        )

        # Final safety rounding by climate variable.
        climate_fields = ["min", "max"]
        if include_mean:
            climate_fields.insert(1, "mean")

        for climate_field in climate_fields:
            if climate_field not in combined_csv.columns:
                continue

            combined_csv[
                climate_field
            ] = pd.to_numeric(
                combined_csv[
                    climate_field
                ],
                errors="coerce",
            )

            for variable_name in (
                combined_csv[
                    "variable"
                ]
                .dropna()
                .unique()
            ):
                variable_mask = (
                    combined_csv[
                        "variable"
                    ]
                    == variable_name
                )

                combined_csv.loc[
                    variable_mask,
                    climate_field,
                ] = (
                    combined_csv.loc[
                        variable_mask,
                        climate_field,
                    ]
                    .round(
                        _output_decimals(
                            variable_name
                        )
                    )
                )

        combined_csv.to_csv(
            output_csv,
            index=False,
        )

        written[
            "csv"
        ] = output_csv

        _log(
            f"[WRITE CSV] {output_csv} | "
            f"records={len(combined_csv)}",
            verbose=verbose,
        )

    _log("=" * 100, verbose=verbose)
    _log(f"[DONE] {output_gpkg}", verbose=verbose)
    _log(
        f"[LAYERS] {', '.join(sorted(written))}",
        verbose=verbose,
    )

    return written