from __future__ import annotations

"""
Flood hazard module for the NARCliM2 workflow package.

This module implements the flood methodology's approximate climate-adjusted flood-frequency
method. It is NOT a hydraulic flood model and does not calculate flood extent,
depth, velocity, drainage-network behaviour, or channel routing.

Hazard metric
-------------
For a selected current design event (for example a 1-in-10-year event), the
module calculates current rainfall excess/runoff, then finds the future event
frequency that produces approximately the same runoff after ARR climate-change
adjustments to rainfall and losses.

Data sources
------------
* Bureau of Meteorology 2016 Design Rainfall / IFD data (current rainfall).
* Australian Rainfall and Runoff (ARR) Data Hub (storm losses and climate
  adjustment factors).

Spatial input
-------------
A point shapefile/GeoPackage. Any number of points is supported. A field must
identify each point as Urban or Rural. The default configuration uses:
* Urban: 1-hour event duration.
* Rural: 24-hour event duration.

ARR standard storm losses are rural losses. Urban losses should ideally be
replaced with locally justified values using the urban loss override arguments.

Folders
-------
Downloaded/raw inputs:
    <output_root>/NARCliM Data/Flood/

Final outputs:
    <output_root>/Climate_Indicies/Hazards/Flood/
"""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import math
import re
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup

FLOOD_CODE_VERSION = "2026-09-10-flood-method-v6"
ARR_DATA_HUB_URL = "https://data.arr-software.org/"
BOM_IFD_URL = "https://www.bom.gov.au/water/designRainfalls/revised-ifd/"

ARR_SCENARIO_NAMES = {
    "ssp126": "SSP1-2.6",
    "ssp245": "SSP2-4.5",
    "ssp370": "SSP3-7.0",
    "ssp585": "SSP5-8.5",
}


@dataclass(frozen=True)
class FloodPoint:
    point_id: str
    name: str
    context: str
    longitude: float
    latitude: float


@dataclass(frozen=True)
class LossParameters:
    initial_loss_mm: float
    continuing_loss_mm_per_hr: float
    source: str


def _log(message: str, *, verbose: bool) -> None:
    if verbose:
        print(message, flush=True)


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")


def _representative_year(start_year: int, end_year: int) -> int:
    """Return midpoint year consistent with the rest of the package."""
    return int(round((int(start_year) + int(end_year)) / 2.0))


def _flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    if isinstance(result.columns, pd.MultiIndex):
        result.columns = [
            " ".join(
                str(part).strip()
                for part in col
                if str(part).strip() and not str(part).startswith("Unnamed:")
            ).strip()
            for col in result.columns
        ]
    else:
        result.columns = [str(c).strip() for c in result.columns]
    return result


def _numeric(value) -> float:
    if value is None:
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    match = re.search(r"-?\d+(?:\.\d+)?", str(value).replace(",", ""))
    return np.nan if match is None else float(match.group(0))


# =============================================================================
# LOCATION INPUT
# =============================================================================

def _infer_context(value: str) -> str:
    text = str(value).strip().lower()

    if "urban" in text:
        return "urban"

    if "rural" in text:
        return "rural"

    raise ValueError(
        f"Could not infer Urban/Rural context from {value!r}. "
        "The context must contain 'Urban' or 'Rural'."
    )


def _read_flood_points(
    *,
    location_path: str | Path,
    name_field: str,
    context_field: str,
) -> tuple[gpd.GeoDataFrame, list[FloodPoint]]:
    """
    Read flood locations from a point shapefile or GeoPackage.
    """

    locations = gpd.read_file(location_path)

    if locations.empty:
        raise ValueError(
            f"No flood locations found: {location_path}"
        )

    if locations.crs is None:
        raise ValueError(
            "Flood location file has no CRS."
        )

    for field in (name_field, context_field):
        if field not in locations.columns:
            raise KeyError(
                f"Field {field!r} not found. "
                f"Available: {list(locations.columns)}"
            )

    if not locations.geometry.geom_type.isin(["Point"]).all():
        raise ValueError(
            "Flood location input must contain Point features only."
        )

    # Convert all locations to latitude/longitude.
    wgs84 = locations.to_crs(4326).copy()

    points: list[FloodPoint] = []

    for n, (_, row) in enumerate(
        wgs84.iterrows(),
        start=1,
    ):

        context = _infer_context(
            row[context_field]
        )

        point_id = (
            f"F{n}_"
            f"{_safe_name(row[name_field])}"
        )

        points.append(
            FloodPoint(
                point_id=point_id,
                name=str(row[name_field]),
                context=context,
                longitude=float(row.geometry.x),
                latitude=float(row.geometry.y),
            )
        )

    wgs84["_flood_point_id"] = [
        p.point_id
        for p in points
    ]

    wgs84["_flood_context"] = [
        p.context
        for p in points
    ]

    return wgs84, points


def _read_coordinate_flood_points(
    *,
    coordinate_points: list[dict],
) -> tuple[gpd.GeoDataFrame, list[FloodPoint]]:
    """
    Create flood locations directly from user-supplied coordinates.

    Each point must contain:

        {
            "name": "Canberra Urban",
            "context": "Urban",
            "latitude": -35.2809,
            "longitude": 149.1300,
        }

    Any number of points can be supplied.

    Notes
    -----
    Coordinates must be supplied in geographic latitude/longitude
    coordinates using WGS84 (EPSG:4326).
    """

    if not coordinate_points:
        raise ValueError(
            "coordinate_points was supplied but contains no points."
        )

    points: list[FloodPoint] = []

    records = []

    # ==================================================================
    # PROCESS USER-DEFINED COORDINATES
    # ==================================================================

    for n, point in enumerate(
        coordinate_points,
        start=1,
    ):

        # --------------------------------------------------------------
        # Check required inputs.
        # --------------------------------------------------------------

        required_fields = {
            "name",
            "context",
            "latitude",
            "longitude",
        }

        missing_fields = (
            required_fields
            - set(point)
        )

        if missing_fields:
            raise ValueError(
                f"Coordinate point {n} is missing "
                f"required field(s): "
                f"{sorted(missing_fields)}"
            )

        # --------------------------------------------------------------
        # Read point information.
        # --------------------------------------------------------------

        name = str(
            point["name"]
        )

        context = _infer_context(
            point["context"]
        )

        latitude = float(
            point["latitude"]
        )

        longitude = float(
            point["longitude"]
        )

        # --------------------------------------------------------------
        # Validate latitude / longitude.
        # --------------------------------------------------------------

        if not -90.0 <= latitude <= 90.0:
            raise ValueError(
                f"Invalid latitude for {name!r}: "
                f"{latitude}. "
                "Latitude must be between -90 and 90."
            )

        if not -180.0 <= longitude <= 180.0:
            raise ValueError(
                f"Invalid longitude for {name!r}: "
                f"{longitude}. "
                "Longitude must be between -180 and 180."
            )

        # --------------------------------------------------------------
        # Create the same FloodPoint structure used by shapefile input.
        #
        # This is important because everything downstream can therefore
        # remain unchanged.
        # --------------------------------------------------------------

        point_id = (
            f"F{n}_"
            f"{_safe_name(name)}"
        )

        flood_point = FloodPoint(
            point_id=point_id,
            name=name,
            context=context,
            longitude=longitude,
            latitude=latitude,
        )

        points.append(
            flood_point
        )

        # --------------------------------------------------------------
        # Create the corresponding spatial record.
        #
        # This lets the final output remain a GeoDataFrame exactly like
        # the shapefile-based workflow.
        # --------------------------------------------------------------

        records.append(
            {
                "Name": name,
                "_flood_point_id": point_id,
                "_flood_context": context,
                "geometry": gpd.points_from_xy(
                    [longitude],
                    [latitude],
                )[0],
            }
        )

    # ==================================================================
    # CREATE WGS84 GEODATAFRAME
    # ==================================================================

    locations = gpd.GeoDataFrame(
        records,
        geometry="geometry",
        crs="EPSG:4326",
    )

    return locations, points


def _load_flood_points(
    *,
    location_path: str | Path | None,
    coordinate_points: list[dict] | None,
    name_field: str,
    context_field: str,
) -> tuple[gpd.GeoDataFrame, list[FloodPoint]]:
    """
    Load flood locations from either user-defined coordinates or a
    spatial file.

    Input priority
    --------------
    1. If coordinate_points is supplied, use those coordinates.
    2. Otherwise, use location_path.

    The returned objects are identical in structure regardless of which
    input method is used, so the flood calculations do not need to know
    where the locations came from.
    """

    # ==================================================================
    # OPTION 1: USER-DEFINED LATITUDE / LONGITUDE
    # ==================================================================

    if coordinate_points is not None:

        if len(coordinate_points) == 0:
            raise ValueError(
                "coordinate_points was supplied but is empty."
            )

        return _read_coordinate_flood_points(
            coordinate_points=coordinate_points,
        )

    # ==================================================================
    # OPTION 2: SHAPEFILE / GEOPACKAGE
    # ==================================================================

    if location_path is not None:

        return _read_flood_points(
            location_path=location_path,
            name_field=name_field,
            context_field=context_field,
        )

    # ==================================================================
    # NO LOCATION INPUT
    # ==================================================================

    raise ValueError(
        "No flood locations were supplied. "
        "Provide either:\n"
        "  1. coordinate_points; or\n"
        "  2. location_path."
    )


# =============================================================================
# DOWNLOAD ARR + BOM INPUTS
# =============================================================================

def _download_arr_html(*, longitude: float, latitude: float, timeout: int) -> tuple[str, str]:
    """Download ARR Data Hub text output using the same request style as the flood methodology."""
    response = requests.get(
        ARR_DATA_HUB_URL,
        params={
            "lon_coord": longitude,
            "lat_coord": latitude,
            "All": 1,
            "type": "text",
        },
        timeout=timeout,
    )
    response.raise_for_status()
    return response.text, response.url


def _download_bom_design_rainfall_html(
    *,
    longitude: float,
    latitude: float,
    design: str,
    timeout: int,
) -> tuple[str, str]:
    """
    Download one BoM Design Rainfall range for a coordinate.

    Parameters
    ----------
    design
        One of:

        "very_frequent"
            Very frequent design rainfall.

        "ifds"
            Standard frequent/infrequent IFD rainfall.

        "rare"
            Rare design rainfall.

    Returns
    -------
    html, url
    """

    valid_designs = {
        "very_frequent",
        "ifds",
        "rare",
    }

    if design not in valid_designs:
        raise ValueError(
            f"Unknown BoM design rainfall range: {design!r}. "
            f"Choose from {sorted(valid_designs)}."
        )

    session = requests.Session()

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/151.0.0.0 Safari/537.36"
        ),
        "Accept": (
            "text/html,application/xhtml+xml,application/xml;q=0.9,"
            "image/avif,image/webp,*/*;q=0.8"
        ),
        "Accept-Language":
            "en-AU,en;q=0.9",

        "Referer":
            BOM_IFD_URL,

        "Connection":
            "keep-alive",
    }

    # Open the base page first to establish the normal web session.
    base_response = session.get(
        BOM_IFD_URL,
        headers=headers,
        timeout=timeout,
    )

    base_response.raise_for_status()

    params = {
        "coordinate_type":
            "dd",

        "design":
            design,

        "latitude":
            latitude,

        "longitude":
            longitude,

        # Request standard minute/hour/day durations.
        "sdmin":
            "true",

        "sdhr":
            "true",

        "sdday":
            "true",

        # Flood runoff calculation requires rainfall depth.
        "values":
            "depths",

        "year":
            "2016",
    }

    response = session.get(
        BOM_IFD_URL,
        params=params,
        headers=headers,
        timeout=timeout,
    )

    response.raise_for_status()

    return (
        response.text,
        response.url,
    )


def _download_all_bom_design_rainfall(
    *,
    longitude: float,
    latitude: float,
    timeout: int,
) -> dict[str, tuple[str, str]]:
    """
    Download all BoM design-rainfall probability ranges.

    Returns
    -------
    dict

    Example
    -------
    {
        "very_frequent": (html, url),
        "ifds": (html, url),
        "rare": (html, url),
    }
    """

    results = {}

    for design in (
        "very_frequent",
        "ifds",
        "rare",
    ):

        html, url = (
            _download_bom_design_rainfall_html(
                longitude=longitude,
                latitude=latitude,
                design=design,
                timeout=timeout,
            )
        )

        results[
            design
        ] = (
            html,
            url,
        )

    return results

def _parse_all_bom_design_rainfall(
    bom_html_by_design: dict[str, tuple[str, str]],
) -> pd.DataFrame:
    """
    Parse and combine Very Frequent, standard IFD, and Rare BoM
    design-rainfall tables.

    Returns
    -------
    DataFrame

    Columns include:
        duration_hours
        event_rate_per_year
        rainfall_depth_mm
        probability_label
        design_range
        source_url
    """

    frames = []

    for (
        design_range,
        (
            html,
            url,
        ),
    ) in bom_html_by_design.items():

        parsed = _parse_bom_ifd(
            html
        ).copy()

        parsed[
            "design_range"
        ] = design_range

        parsed[
            "source_url"
        ] = url

        frames.append(
            parsed
        )

    if not frames:
        raise ValueError(
            "No BoM design-rainfall tables were parsed."
        )

    combined = pd.concat(
        frames,
        ignore_index=True,
    )

    # Remove duplicates around range boundaries.
    #
    # For example 1 EY / 1% may appear in adjacent products.
    combined = (
        combined
        .sort_values(
            [
                "duration_hours",
                "event_rate_per_year",
            ]
        )
        .drop_duplicates(
            subset=[
                "duration_hours",
                "event_rate_per_year",
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    return combined

# =============================================================================
# ARR PARSERS
# =============================================================================

def _parse_arr_storm_losses(raw: str) -> LossParameters:
    """Parse ARR storm initial and continuing losses from saved ARR input.

    The workflow may reuse either the older HTML response or flood-method
    plain-text ARR output. This parser supports both formats without changing
    the existing ARR download/reuse workflow.
    """
    text = _arr_plain_text(raw)

    def _find_value(label: str) -> float | None:
        # flood-method ARR text is comma separated, for example:
        # Storm Initial Losses (mm), 25.0
        for line in text.splitlines():
            if label.lower() in line.lower():
                parts = line.split(",")
                if len(parts) >= 2:
                    value = _numeric(parts[1])
                    if np.isfinite(value):
                        return float(value)

        # Backwards-compatible fallback for older HTML/plain-text layouts.
        label_pattern = re.escape(label).replace(r"\ ", r"\s+")
        match = re.search(
            label_pattern + r"\s*[,=:]?\s*(-?\d+(?:\.\d+)?)",
            text,
            flags=re.I,
        )
        if match is not None:
            return float(match.group(1))
        return None

    il = _find_value("Storm Initial Losses (mm)")
    cl = _find_value("Storm Continuing Losses (mm/h)")

    if il is None or cl is None:
        raise ValueError(
            "Could not parse ARR storm initial/continuing losses from the saved ARR input."
        )

    return LossParameters(
        initial_loss_mm=float(il),
        continuing_loss_mm_per_hr=float(cl),
        source="ARR Data Hub storm losses",
    )


def _table_df(table_tag) -> pd.DataFrame:
    return _flatten_columns(pd.read_html(StringIO(str(table_tag)))[0])


def _scenario_from_heading(text: str) -> str | None:
    compact = re.sub(r"\s+", "", str(text)).upper()
    for package_name, arr_name in ARR_SCENARIO_NAMES.items():
        if re.sub(r"\s+", "", arr_name).upper() in compact:
            return package_name
    return None


def _parse_arr_rainfall_factor_tables(html: str) -> dict[str, pd.DataFrame]:
    """Parse ARR rainfall factors by SSP."""
    soup = BeautifulSoup(html, "html.parser")
    output: dict[str, pd.DataFrame] = {}
    for table in soup.find_all("table"):
        heading = table.find_previous(["h5", "h6"])
        if heading is None:
            continue
        scenario = _scenario_from_heading(heading.get_text(" ", strip=True))
        if scenario is None:
            continue
        # Current ARR page structure places SSP rainfall tables under the
        # nearest preceding h4 heading "Rainfall Factors".
        section = heading.find_previous("h4")
        if section is None or "rainfall" not in section.get_text(" ", strip=True).lower():
            continue
        df = _table_df(table)
        year_col = next((c for c in df.columns if "year" in c.lower()), None)
        if year_col is None:
            continue
        df = df.rename(columns={year_col: "Year"})
        df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
        output[scenario] = df.dropna(subset=["Year"])
    if not output:
        raise ValueError("Could not parse ARR rainfall climate-change factor tables.")
    return output


def _parse_arr_loss_factor_tables(html: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    soup = BeautifulSoup(html, "html.parser")
    initial = continuing = None
    for table in soup.find_all("table"):
        heading = table.find_previous(["h5", "h6"])
        if heading is None:
            continue
        text = heading.get_text(" ", strip=True).lower()
        if "initial loss" in text and "adjustment" in text:
            initial = _table_df(table)
        elif "continuing loss" in text and "adjustment" in text:
            continuing = _table_df(table)
    if initial is None or continuing is None:
        raise ValueError("Could not parse ARR loss adjustment tables.")
    return initial, continuing


def _year_column(df: pd.DataFrame) -> str:
    return next((c for c in df.columns if "year" in str(c).lower()), df.columns[0])


def _scenario_column(df: pd.DataFrame, scenario: str) -> str:
    target = re.sub(r"\s+", "", ARR_SCENARIO_NAMES[scenario]).lower()
    for col in df.columns:
        if target in re.sub(r"\s+", "", str(col)).lower():
            return col
    raise KeyError(f"ARR loss-factor column not found for {scenario}. Columns={list(df.columns)}")


def _interpolate_table_year(*, df: pd.DataFrame, value_column: str, year: int) -> float:
    yc = _year_column(df)
    years = pd.to_numeric(df[yc], errors="coerce").to_numpy(float)
    values = pd.to_numeric(df[value_column], errors="coerce").to_numpy(float)
    valid = np.isfinite(years) & np.isfinite(values)
    years, values = years[valid], values[valid]
    order = np.argsort(years)
    years, values = years[order], values[order]
    if not len(years) or year < years.min() or year > years.max():
        raise ValueError(f"Requested climate year {year} is outside the ARR factor table.")
    return float(np.interp(year, years, values))


def _duration_anchor_hours(label: str) -> float | None:
    text = str(label).strip().lower()
    if "<1" in text and "hour" in text:
        return 1.0
    if ">24" in text and "hour" in text:
        return 24.0
    if "hour" not in text:
        return None
    m = re.search(r"(\d+(?:\.\d+)?)", text)
    return None if m is None else float(m.group(1))


def _rainfall_factor_for_duration(
    *, df: pd.DataFrame, representative_year: int, duration_hours: float
) -> float:
    yc = next((c for c in df.columns if "year" in str(c).lower()), None)
    if yc is None:
        raise ValueError("ARR rainfall-factor table has no Year column.")
    pairs = []
    for col in df.columns:
        if col == yc:
            continue
        anchor = _duration_anchor_hours(col)
        if anchor is not None:
            pairs.append((anchor, col))
    pairs.sort()
    if not pairs:
        raise ValueError("No ARR rainfall-factor duration columns found.")
    hours = np.array([p[0] for p in pairs], dtype=float)
    factors = np.array(
        [_interpolate_table_year(df=df, value_column=p[1], year=representative_year) for p in pairs],
        dtype=float,
    )
    # the flood methodology's current simplified implementation uses 1h and 24h, which map
    # to the endpoint classes. Intermediate durations are linearly interpolated.
    d = float(np.clip(duration_hours, hours.min(), hours.max()))
    return float(np.interp(d, hours, factors))



# =============================================================================
# FLOOD ARR TEXT PARSERS / CLIMATE-CHANGE FACTORS
# =============================================================================

def _arr_plain_text(raw: str) -> str:
    """Return ARR response as plain text whether the response is text or HTML."""
    if "<html" in raw.lower() or "<body" in raw.lower():
        return BeautifulSoup(raw, "html.parser").get_text("\n", strip=False)
    return raw


def _find_arr_section(raw: str, start_marker: str, end_marker: str) -> list[str]:
    """Extract lines between ARR text-file section markers, matching the flood methodology's helper."""
    lines: list[str] = []
    reading = False
    for line in _arr_plain_text(raw).splitlines():
        if end_marker in line:
            break
        if reading:
            lines.append(line)
        if start_marker in line:
            reading = True
            continue
    if not lines:
        raise ValueError(
            f"ARR section {start_marker!r} to {end_marker!r} was not found. "
            "Use overwrite=True once if an older HTML ARR cache is being reused."
        )
    return lines


def _parse_arr_text_table(raw: str, start_marker: str, end_marker: str) -> pd.DataFrame:
    """Parse one comma-separated ARR text section into a DataFrame."""
    lines = _find_arr_section(raw, start_marker, end_marker)
    df = pd.read_csv(StringIO("\n".join(lines)), header=0)
    if df.empty:
        raise ValueError(f"ARR section {start_marker!r} parsed to an empty table.")
    return _flatten_columns(df)


def _parse_arr_temperature_changes(raw: str) -> pd.DataFrame:
    return _parse_arr_text_table(
        raw,
        "[TEMPERATURE_CHANGES]",
        "[END_TEMPERATURE_CHANGES]",
    )


def _parse_arr_loss_tables(raw: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    initial = _parse_arr_text_table(
        raw,
        "[Climate_Change_INITIAL_LOSS]",
        "[END_Climate_Change_INITIAL_LOSS]",
    )
    continuing = _parse_arr_text_table(
        raw,
        "[Climate_Change_CONTINUING_LOSS]",
        "[END_Climate_Change_CONTINUING_LOSS]",
    )
    return initial, continuing


def _table_year_values(df: pd.DataFrame) -> tuple[str, np.ndarray]:
    """Return the most likely year column and numeric year values."""
    year_col = next(
        (c for c in df.columns if "year" in str(c).lower()),
        df.columns[0],
    )
    years = pd.to_numeric(df[year_col], errors="coerce").to_numpy(float)
    return year_col, years


def _arr_table_value(df: pd.DataFrame, *, year: int, column: str) -> float:
    """Read/interpolate a the flood methodology ARR table value for the requested horizon year."""
    if column not in df.columns:
        # Be forgiving about whitespace in ARR headings.
        compact_target = re.sub(r"\s+", "", column).lower()
        match = next(
            (
                c for c in df.columns
                if compact_target in re.sub(r"\s+", "", str(c)).lower()
            ),
            None,
        )
        if match is None:
            raise KeyError(f"ARR column {column!r} not found. Columns={list(df.columns)}")
        column = match

    _, years = _table_year_values(df)
    values = pd.to_numeric(df[column], errors="coerce").to_numpy(float)
    valid = np.isfinite(years) & np.isfinite(values)
    years = years[valid]
    values = values[valid]
    if years.size == 0:
        raise ValueError(f"No numeric values found in ARR column {column!r}.")

    order = np.argsort(years)
    years = years[order]
    values = values[order]

    if year < years.min() or year > years.max():
        raise ValueError(
            f"Requested year {year} is outside ARR table range "
            f"{years.min():g}-{years.max():g}."
        )
    return float(np.interp(year, years, values))


def _arr_delta_t(
    temperature_table: pd.DataFrame,
    *,
    year: int,
    scenario: str,
) -> float:
    """Temperature increase used by the flood methodology for the horizon/scenario."""
    arr_scenario = ARR_SCENARIO_NAMES[scenario]
    return _arr_table_value(
        temperature_table,
        year=year,
        column=arr_scenario,
    )


def _arr_loss_per_degree_factor(
    loss_table: pd.DataFrame,
    *,
    year: int,
    scenario: str,
) -> float:
    """ARR per-degree loss factor used in the flood methodology's ``factor ** deltaT`` formula."""
    arr_scenario = ARR_SCENARIO_NAMES[scenario]
    return _arr_table_value(
        loss_table,
        year=year,
        column=f"Losses {arr_scenario}",
    )


def _rainfall_alpha_for_duration(duration_hours: float) -> float:
    """Rainfall percentage increase per degree C used in the flood methodology's calculation.

    Exact method defaults:
        1 hour  -> 15 % / degree C
        24 hour ->  8 % / degree C

    The workflow defaults remain 1 h (urban) and 24 h (rural), so those are
    exact reproductions. For any optional intermediate duration supplied by a
    user, alpha is linearly interpolated only to keep the function usable.
    """
    d = float(duration_hours)
    if math.isclose(d, 1.0, abs_tol=1e-9):
        return 15.0
    if math.isclose(d, 24.0, abs_tol=1e-9):
        return 8.0
    return float(np.interp(np.clip(d, 1.0, 24.0), [1.0, 24.0], [15.0, 8.0]))


def _rainfall_climate_factor(*, duration_hours: float, delta_t: float) -> float:
    """the flood methodology's climate-adjusted rainfall multiplier: (1 + alpha/100) ** deltaT."""
    alpha = _rainfall_alpha_for_duration(duration_hours)
    return float((1.0 + alpha / 100.0) ** float(delta_t))


def _reference_event_rate(current_ari_year: float) -> float:
    """Return the flood methodology's reference EY value for 1-, 10-, and 100-year events."""
    ari = float(current_ari_year)
    if math.isclose(ari, 1.0, abs_tol=1e-9):
        return 1.0
    if math.isclose(ari, 10.0, abs_tol=1e-9):
        return 0.11
    if math.isclose(ari, 100.0, abs_tol=1e-9):
        return 0.01

    warnings.warn(
        f"ARI {ari:g} years is not one of the flood methodology's 1/10/100-year reference events; "
        "using 1/ARI for this optional custom event.",
        RuntimeWarning,
    )
    return 1.0 / ari


# =============================================================================
# BOM IFD PARSER
# =============================================================================

def _duration_to_hours(value) -> float | None:
    text = str(value).lower().strip()
    m = re.search(r"(\d+(?:\.\d+)?)", text)
    if m is None:
        return None
    n = float(m.group(1))
    if "min" in text:
        return n / 60.0
    if "hour" in text or "hr" in text:
        return n
    if "day" in text:
        return n * 24.0
    return None


def _event_rate_from_label(label: str) -> float | None:
    """Convert BoM EY/AEP/ARI column label to expected events per year."""
    text = re.sub(r"\s+", " ", str(label)).strip().lower()
    m = re.search(r"(\d+(?:\.\d+)?)\s*ey", text)
    if m:
        return float(m.group(1))
    m = re.search(r"1\s*(?:in|/)\s*(\d+(?:\.\d+)?)", text)
    if m:
        return 1.0 / float(m.group(1))
    m = re.search(r"(\d+(?:\.\d+)?)\s*%", text)
    if m:
        p = float(m.group(1)) / 100.0
        if 0 < p < 1:
            # Poisson conversion; 63.2% AEP ~= 1 EY.
            return -math.log(1.0 - p)
    return None


def _prepare_html_table_for_pandas(
    table_tag,
) -> str:
    """
    Clean one HTML table before passing it to pandas.read_html().

    Why this is needed
    ------------------
    The BoM Design Rainfall webpage contains some HTML layout tables with
    non-standard colspan/rowspan attributes such as:

        colspan="100%"

    pandas expects colspan and rowspan to contain integer values only.
    Therefore attempting to parse the complete BoM page directly can raise:

        ValueError:
        invalid literal for int() with base 10: '100%'

    This helper removes invalid colspan/rowspan values while retaining
    legitimate integer spans.

    The cleaning is applied only to an in-memory copy used for parsing;
    the original downloaded BoM HTML file remains unchanged.
    """

    # Convert to text and reparse so that we do not modify the original
    # BeautifulSoup object.
    table_soup = BeautifulSoup(
        str(table_tag),
        "html.parser",
    )

    for element in table_soup.find_all(
        True
    ):

        # --------------------------------------------------------------
        # Clean colspan.
        # --------------------------------------------------------------

        if "colspan" in element.attrs:

            colspan = str(
                element.attrs[
                    "colspan"
                ]
            ).strip()

            # Valid HTML table spans for pandas must be positive integers.
            if not colspan.isdigit():

                # Example:
                #     colspan="100%"
                #
                # This is generally webpage layout formatting rather than
                # meaningful tabular structure, so remove the attribute.
                del element.attrs[
                    "colspan"
                ]

        # --------------------------------------------------------------
        # Clean rowspan.
        # --------------------------------------------------------------

        if "rowspan" in element.attrs:

            rowspan = str(
                element.attrs[
                    "rowspan"
                ]
            ).strip()

            if not rowspan.isdigit():

                del element.attrs[
                    "rowspan"
                ]

    return str(
        table_soup
    )

def _event_rate_from_probability_label(
    label: str,
) -> float | None:
    """Convert BoM probability/frequency labels using the flood methodology's EY convention.

    the flood methodology's calculation files use the following practical mapping when the
    normal IFD table is merged with the very-frequent and rare products:

        63.2% -> 1 EY
        50%   -> 0.5 EY
        20%   -> 0.2 EY
        10%   -> 0.11 EY
        5%    -> 0.05 EY
        2%    -> 0.02 EY
        1%    -> 0.01 EY

    Rare-event labels such as ``1 in 200`` remain 1 / ARI, while columns
    already labelled in EY are read directly.
    """
    if label is None:
        return None

    text = re.sub(r"\s+", " ", str(label).strip().lower())

    # 1) Explicit EY columns: 12EY, 6EY, 0.5EY, etc.
    ey_match = re.search(r"(\d+(?:\.\d+)?)\s*ey", text)
    if ey_match is not None:
        return float(ey_match.group(1))

    # 2) Rare-event ARI labels: 1 in 200, 1/500, etc.
    ari_match = re.search(r"1\s*(?:in|/)\s*(\d+(?:\.\d+)?)", text)
    if ari_match is not None:
        ari_years = float(ari_match.group(1))
        return None if ari_years <= 0 else 1.0 / ari_years

    # 3) the flood methodology's mapping from AEP columns to EY column names.
    percent_match = re.search(r"(\d+(?:\.\d+)?)\s*%", text)
    if percent_match is not None:
        p = float(percent_match.group(1))
        ey_map = {
            63.2: 1.0,
            50.0: 0.5,
            20.0: 0.2,
            10.0: 0.11,
            5.0: 0.05,
            2.0: 0.02,
            1.0: 0.01,
        }
        for key, value in ey_map.items():
            if math.isclose(p, key, rel_tol=0.0, abs_tol=1e-6):
                return value

    return None


def _parse_bom_ifd(
    html: str,
) -> pd.DataFrame:
    """
    Parse BoM 2016 IFD design-rainfall depths into long format.

    Parameters
    ----------
    html
        HTML returned from the BoM Design Rainfall webpage.

    Returns
    -------
    pandas.DataFrame
        Long-format table containing:

        - duration_hours
        - event_rate_per_year
        - rainfall_depth_mm
        - probability_label

    Processing approach
    -------------------
    The BoM webpage contains several HTML tables, including webpage-layout
    tables that are unrelated to design rainfall.

    Rather than passing the complete page directly to ``pd.read_html()``,
    this function:

        1. finds HTML tables individually;
        2. cleans malformed colspan/rowspan attributes;
        3. attempts to parse each table separately;
        4. skips tables that are not valid tabular data;
        5. identifies tables containing rainfall duration and probability
           information;
        6. converts the BoM table to a consistent long-format dataset.

    This prevents unrelated HTML layout tables from causing the entire flood
    workflow to fail.
    """

    soup = BeautifulSoup(
        html,
        "html.parser",
    )

    html_tables = soup.find_all(
        "table"
    )

    if not html_tables:

        raise ValueError(
            "No HTML tables were found in the BoM IFD response."
        )

    records = []

    parsed_table_count = 0
    skipped_table_count = 0

    # ==================================================================
    # PROCESS EACH HTML TABLE SEPARATELY
    # ==================================================================

    for table_number, table_tag in enumerate(
        html_tables,
        start=1,
    ):

        try:

            cleaned_table_html = (
                _prepare_html_table_for_pandas(
                    table_tag
                )
            )

            parsed_tables = pd.read_html(
                StringIO(
                    cleaned_table_html
                ),

                # lxml is generally robust for BoM HTML.
                flavor="lxml",

                # Some useful result tables can be hidden/shown by the
                # webpage interface. We still want their underlying data.
                displayed_only=False,
            )

        except Exception:

            # ----------------------------------------------------------
            # Some BoM webpage tables are purely for page layout.
            # They are irrelevant to the IFD data and should not stop
            # processing.
            # ----------------------------------------------------------

            skipped_table_count += 1

            continue

        if not parsed_tables:

            skipped_table_count += 1

            continue

        # One HTML <table> normally produces one DataFrame, but iterate
        # safely in case pandas returns more than one.
        for raw_table in parsed_tables:

            parsed_table_count += 1

            df = _flatten_columns(
                raw_table
            )

            if df.empty:

                continue

            # ==========================================================
            # IDENTIFY THE DURATION COLUMN
            #
            # BoM duration cells can look like:
            #
            #     60 min
            #     1 hour
            #     2 hour
            #     24 hour
            #     1 day
            #
            # We inspect the values rather than relying on an exact column
            # heading because BoM HTML formatting can vary.
            # ==========================================================

            duration_column = None

            for column in df.columns:

                sample_values = (
                    df[
                        column
                    ]
                    .dropna()
                    .astype(
                        str
                    )
                    .head(
                        20
                    )
                    .tolist()
                )

                if not sample_values:

                    continue

                parsed_durations = [
                    _duration_to_hours(
                        value
                    )
                    for value in sample_values
                ]

                valid_duration_count = sum(
                    duration is not None
                    for duration
                    in parsed_durations
                )

                # Require more than one recognisable duration so that an
                # unrelated webpage table is not accidentally interpreted
                # as an IFD table.
                if valid_duration_count >= 2:

                    duration_column = column

                    break

            if duration_column is None:

                continue

            # ==========================================================
            # IDENTIFY PROBABILITY / FREQUENCY COLUMNS
            #
            # Examples include:
            #
            #     12 EY
            #     2 EY
            #     1 EY
            #     50%
            #     20%
            #     10%
            #     5%
            #     2%
            #     1%
            #     1 in 100
            #
            # _event_rate_from_probability_label() converts all of these
            # to a common expected-events-per-year representation.
            # ==========================================================

            probability_columns = []

            for column in df.columns:

                if column == duration_column:

                    continue

                event_rate = (
                    _event_rate_from_probability_label(
                        column
                    )
                )

                if event_rate is None:

                    continue

                probability_columns.append(
                    (
                        column,
                        event_rate,
                    )
                )

            if not probability_columns:

                continue

            # ==========================================================
            # CONVERT TABLE FROM WIDE → LONG FORMAT
            # ==========================================================

            for _, row in df.iterrows():

                duration_hours = (
                    _duration_to_hours(
                        row[
                            duration_column
                        ]
                    )
                )

                if duration_hours is None:

                    continue

                for (
                    probability_column,
                    event_rate,
                ) in probability_columns:

                    rainfall_depth = (
                        _numeric(
                            row[
                                probability_column
                            ]
                        )
                    )

                    if not np.isfinite(
                        rainfall_depth
                    ):

                        continue

                    # Rainfall depth should never be negative.
                    if rainfall_depth < 0:

                        continue

                    records.append(
                        {
                            "duration_hours":
                                float(
                                    duration_hours
                                ),

                            "event_rate_per_year":
                                float(
                                    event_rate
                                ),

                            "rainfall_depth_mm":
                                float(
                                    rainfall_depth
                                ),

                            "probability_label":
                                str(
                                    probability_column
                                ),

                            # Useful during QA to know which HTML table
                            # supplied the value.
                            "source_table":
                                table_number,
                        }
                    )

    # ==================================================================
    # VALIDATE RESULTS
    # ==================================================================

    if not records:

        raise ValueError(
            "BoM webpage was downloaded successfully, but no IFD "
            "design-rainfall values could be parsed. "
            f"HTML tables found={len(html_tables)}, "
            f"tables parsed={parsed_table_count}, "
            f"tables skipped={skipped_table_count}."
        )

    result = pd.DataFrame(
        records
    )

    # ==================================================================
    # REMOVE DUPLICATES
    #
    # BoM can contain repeated representations of the same result within
    # the webpage. Keep a single depth for each duration/frequency pair.
    # ==================================================================

    result = (
        result
        .sort_values(
            [
                "duration_hours",
                "event_rate_per_year",
                "source_table",
            ]
        )
        .drop_duplicates(
            subset=[
                "duration_hours",
                "event_rate_per_year",
            ],
            keep="first",
        )
        .sort_values(
            [
                "duration_hours",
                "event_rate_per_year",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    return result


def _ifd_for_duration(ifd: pd.DataFrame, duration_hours: float) -> pd.DataFrame:
    available = np.array(sorted(ifd["duration_hours"].unique()), dtype=float)
    match = available[np.isclose(available, duration_hours, atol=1e-6, rtol=0)]
    if match.size == 0:
        raise ValueError(
            f"Duration {duration_hours:g}h was not found in BoM IFD data. "
            f"Available durations: {available.tolist()}"
        )
    return ifd[np.isclose(ifd["duration_hours"], float(match[0]))].sort_values(
        "event_rate_per_year"
    ).copy()


def _interpolate_ifd_depth(duration_ifd: pd.DataFrame, target_event_rate: float) -> float:
    rates = duration_ifd["event_rate_per_year"].to_numpy(float)
    depths = duration_ifd["rainfall_depth_mm"].to_numpy(float)
    valid = np.isfinite(rates) & np.isfinite(depths) & (rates > 0)
    rates, depths = rates[valid], depths[valid]
    order = np.argsort(rates)
    rates, depths = rates[order], depths[order]
    if target_event_rate < rates.min() or target_event_rate > rates.max():
        raise ValueError(
            f"Target event rate {target_event_rate:.6f}/yr is outside the BoM IFD range "
            f"{rates.min():.6f}-{rates.max():.6f}/yr."
        )
    return float(np.interp(np.log(target_event_rate), np.log(rates), depths))


# =============================================================================
# RAINFALL EXCESS / RUNOFF
# =============================================================================

def calculate_rainfall_excess(
    *,
    rainfall_depth_mm: float,
    duration_hours: float,
    initial_loss_mm: float,
    continuing_loss_mm_per_hr: float,
    timesteps: int = 60,
) -> pd.DataFrame:
    """flood-method compatible multi-timestep IL/CL rainfall-excess calculation.

    This reproduces the logic in the flood methodology's ``basic_hydrograph`` routine:

    1. distribute total rainfall uniformly over the storm timesteps;
    2. convert continuing loss from mm/h to a fixed mm/timestep capacity;
    3. fill the initial-loss store first;
    4. once the initial-loss store is full, apply the FULL timestep continuing
       loss capacity (including the timestep in which IL becomes full);
    5. remaining rainfall is rainfall excess/runoff.

    Note that this deliberately differs from the previous version of this
    module, which prorated continuing loss when IL was exhausted part-way
    through a timestep. the flood methodology's code does not prorate that transition step.
    """
    rainfall_depth_mm = float(rainfall_depth_mm)
    duration_hours = float(duration_hours)
    initial_loss_mm = max(0.0, float(initial_loss_mm))
    continuing_loss_mm_per_hr = max(0.0, float(continuing_loss_mm_per_hr))
    timesteps = int(timesteps)

    if rainfall_depth_mm < 0 or duration_hours <= 0 or timesteps < 1:
        raise ValueError("Invalid rainfall depth, duration, or timestep count.")

    dt = duration_hours / timesteps
    rain_per_step = rainfall_depth_mm / timesteps
    cl_capacity_per_step = continuing_loss_mm_per_hr * dt

    il_full = initial_loss_mm <= 0.0
    cumulative_il = 0.0
    rows = []

    for step in range(1, timesteps + 1):
        previous_il = cumulative_il

        # the flood methodology: il_step is the cumulative IL after adding this timestep's rain,
        # capped at the specified storm initial loss.
        cumulative_il_candidate = min(cumulative_il + rain_per_step, initial_loss_mm)

        # Rainfall remaining after satisfying any outstanding initial loss.
        excess_available_for_cl = max(
            0.0,
            cumulative_il + rain_per_step - initial_loss_mm,
        )

        if cumulative_il_candidate == initial_loss_mm:
            il_full = True
            cumulative_il = initial_loss_mm
        else:
            cumulative_il = cumulative_il + rain_per_step

        # Convert the flood methodology's cumulative IL state into an actual loss for this step
        # so the QA table has a proper water-balance column.
        initial_loss_step = max(0.0, cumulative_il - previous_il)

        if il_full:
            continuing_loss_step = min(
                excess_available_for_cl,
                cl_capacity_per_step,
            )
            runoff_step = excess_available_for_cl - continuing_loss_step
        else:
            continuing_loss_step = 0.0
            runoff_step = 0.0

        rows.append(
            {
                "timestep": step,
                "time_start_hr": (step - 1) * dt,
                "time_end_hr": step * dt,
                "rainfall_mm": rain_per_step,
                "initial_loss_mm": initial_loss_step,
                "continuing_loss_mm": continuing_loss_step,
                "runoff_mm": runoff_step,
                "remaining_initial_loss_mm": max(
                    0.0,
                    initial_loss_mm - cumulative_il,
                ),
                "continuing_loss_capacity_mm": cl_capacity_per_step,
                "initial_loss_full": il_full,
            }
        )

    hydro = pd.DataFrame(rows)

    # Strong QA check: rainfall must equal IL + CL + runoff.
    balance_error = (
        hydro["rainfall_mm"].sum()
        - hydro["initial_loss_mm"].sum()
        - hydro["continuing_loss_mm"].sum()
        - hydro["runoff_mm"].sum()
    )
    if not math.isclose(balance_error, 0.0, abs_tol=1e-8):
        raise RuntimeError(
            "Rainfall-excess water balance failed: "
            f"error={balance_error:.12g} mm"
        )

    return hydro


def _total_runoff(**kwargs) -> float:
    return float(calculate_rainfall_excess(**kwargs)["runoff_mm"].sum())

RUNOFF_MATCH_THRESHOLD_MM = 0.1


def _minimum_rainfall_depth_for_runoff(
    *,
    duration_hours: float,
    initial_loss_mm: float,
    continuing_loss_mm_per_hr: float,
    timesteps: int,
    target_runoff_mm: float = RUNOFF_MATCH_THRESHOLD_MM,
) -> float:
    """Return the minimum uniform storm depth that produces target runoff.

    This uses the same rainfall-excess calculation as the main workflow, so no
    separate hydrologic assumptions are introduced. A bisection search is used
    because total runoff is monotonic with increasing rainfall depth.
    """
    if target_runoff_mm <= 0:
        raise ValueError("target_runoff_mm must be > 0.")

    lower = 0.0
    upper = max(
        1.0,
        float(initial_loss_mm)
        + float(continuing_loss_mm_per_hr) * float(duration_hours)
        + float(target_runoff_mm),
    )

    # Expand the upper bound until the target runoff is reached.
    for _ in range(60):
        if _total_runoff(
            rainfall_depth_mm=upper,
            duration_hours=duration_hours,
            initial_loss_mm=initial_loss_mm,
            continuing_loss_mm_per_hr=continuing_loss_mm_per_hr,
            timesteps=timesteps,
        ) >= target_runoff_mm:
            break
        upper *= 2.0
    else:
        raise RuntimeError(
            "Could not bracket the rainfall depth needed to produce the "
            f"target runoff ({target_runoff_mm:g} mm)."
        )

    # Bisection provides a stable threshold without changing the hydrograph logic.
    for _ in range(80):
        midpoint = 0.5 * (lower + upper)
        runoff = _total_runoff(
            rainfall_depth_mm=midpoint,
            duration_hours=duration_hours,
            initial_loss_mm=initial_loss_mm,
            continuing_loss_mm_per_hr=continuing_loss_mm_per_hr,
            timesteps=timesteps,
        )
        if runoff >= target_runoff_mm:
            upper = midpoint
        else:
            lower = midpoint

    return float(upper)


def _rainfall_shortfall_to_runoff(
    *,
    rainfall_depth_mm: float,
    duration_hours: float,
    initial_loss_mm: float,
    continuing_loss_mm_per_hr: float,
    timesteps: int,
    target_runoff_mm: float = RUNOFF_MATCH_THRESHOLD_MM,
) -> tuple[float, float]:
    """Return (shortfall_mm, threshold_depth_mm) for near-zero-runoff events.

    ``shortfall_mm`` is the additional total storm depth required to reach the
    target runoff.  It is always >= 0 and is therefore safe for the log-log
    """
    threshold_depth = _minimum_rainfall_depth_for_runoff(
        duration_hours=duration_hours,
        initial_loss_mm=initial_loss_mm,
        continuing_loss_mm_per_hr=continuing_loss_mm_per_hr,
        timesteps=timesteps,
        target_runoff_mm=target_runoff_mm,
    )
    shortfall = max(0.0, float(threshold_depth) - float(rainfall_depth_mm))
    return float(shortfall), float(threshold_depth)


def _match_future_event_rate(
    *,
    current_runoff_mm: float,
    current_rainfall_shortfall_mm: float,
    future_curve: pd.DataFrame,
) -> tuple[float, str]:
    """Match current event to a future EY.

    1. If current runoff >= 0.1 mm, use log-log interpolation of
       rainfall excess/runoff against EY.
    2. If current runoff < 0.1 mm, use a corrected rainfall-shortfall metric:
       the extra storm depth required to produce 0.1 mm runoff.  This preserves
       intended fallback while avoiding the NaN/non-positive values in
       the supplied ``r_short`` implementation.
    """
    if not np.isfinite(current_runoff_mm):
        return np.nan, "invalid_current_runoff"

    # ------------------------------------------------------------------
    # NORMAL RUNOFF BRANCH
    # ------------------------------------------------------------------
    if current_runoff_mm >= RUNOFF_MATCH_THRESHOLD_MM:
        curve = future_curve[["event_rate_per_year", "future_runoff_mm"]].dropna().copy()
        curve = curve[
            (curve["event_rate_per_year"] > 0)
            & (curve["future_runoff_mm"] > 0)
        ]
        curve = curve.sort_values("future_runoff_mm").drop_duplicates("future_runoff_mm")

        if len(curve) < 2:
            return np.nan, "insufficient_positive_future_runoff_curve"

        runoff = curve["future_runoff_mm"].to_numpy(float)
        rates = curve["event_rate_per_year"].to_numpy(float)

        if current_runoff_mm < runoff.min() or current_runoff_mm > runoff.max():
            return np.nan, "current_runoff_outside_future_curve"

        log_rate = np.interp(
            np.log10(current_runoff_mm),
            np.log10(runoff),
            np.log10(rates),
        )
        return float(10.0 ** log_rate), "matched_loglog_runoff"

    # ------------------------------------------------------------------
    # NEAR-ZERO RUNOFF / RAINFALL-SHORTFALL BRANCH
    # ------------------------------------------------------------------
    if (
        not np.isfinite(current_rainfall_shortfall_mm)
        or current_rainfall_shortfall_mm <= 0
    ):
        return np.nan, "invalid_current_rainfall_shortfall"

    required = {
        "event_rate_per_year",
        "future_rainfall_shortfall_mm",
    }
    if not required.issubset(future_curve.columns):
        return np.nan, "missing_future_rainfall_shortfall"

    curve = future_curve[
        ["event_rate_per_year", "future_rainfall_shortfall_mm"]
    ].dropna().copy()
    curve = curve[
        (curve["event_rate_per_year"] > 0)
        & (curve["future_rainfall_shortfall_mm"] > 0)
    ]
    curve = (
        curve
        .sort_values("future_rainfall_shortfall_mm")
        .drop_duplicates("future_rainfall_shortfall_mm")
    )

    if len(curve) < 2:
        return np.nan, "insufficient_positive_future_shortfall_curve"

    shortfalls = curve["future_rainfall_shortfall_mm"].to_numpy(float)
    rates = curve["event_rate_per_year"].to_numpy(float)

    if (
        current_rainfall_shortfall_mm < shortfalls.min()
        or current_rainfall_shortfall_mm > shortfalls.max()
    ):
        return np.nan, "current_shortfall_outside_future_curve"

    log_rate = np.interp(
        np.log10(current_rainfall_shortfall_mm),
        np.log10(shortfalls),
        np.log10(rates),
    )
    return float(10.0 ** log_rate), "matched_loglog_shortfall"


# =============================================================================
# URBAN/RURAL LOSS SELECTION
# =============================================================================

def _current_losses_for_context(
    *,
    context: str,
    arr_losses: LossParameters,
    urban_initial_loss_mm: float | None,
    urban_continuing_loss_mm_per_hr: float | None,
) -> tuple[LossParameters, str]:
    """Select losses using the flood methodology's urban/rural assumptions."""
    if context == "rural":
        return arr_losses, "arr_rural_losses"

    # the flood methodology uses fixed urban losses for the 1-hour urban event in both the
    # current and future calculations: IL = 2 mm and CL = 0 mm/h.
    urban_il = 2.0 if urban_initial_loss_mm is None else float(urban_initial_loss_mm)
    urban_cl = (
        0.0
        if urban_continuing_loss_mm_per_hr is None
        else float(urban_continuing_loss_mm_per_hr)
    )
    return (
        LossParameters(
            urban_il,
            urban_cl,
            "urban losses (default IL=2 mm, CL=0 mm/h)",
        ),
        "urban_losses",
    )


# =============================================================================
# FLOOD HAZARD VISUALISATION
# =============================================================================
def plot_flood_frequency_comparison(
    flood_results,
    output_dir,
    time_windows,
    x_max=None,
    y_max=None,
    axis_padding=0.05,
    dpi=300,
    show=False,
):
    """
    Create current-versus-future flood-frequency plots.

    One plot is produced for each scenario / time-horizon combination.

    Main features
    -------------
    1. Automatically detects the current ARIs available in the results.
       For example, if the workflow uses:

           current_aris_years = [5, 10, 100]

       the plotting function automatically uses those events.

    2. Uses one common x-axis maximum and one common y-axis maximum
       across ALL plots.

       This makes the different scenario / horizon figures directly
       comparable.

    3. Allows x_max and y_max to be supplied manually.

    4. Plots each location/context separately rather than averaging
       different locations together.

    5. Adds a 1:1 line:
           above line = more frequent in future
           below line = less frequent in future

    6. Labels each point with its current ARI.

    Parameters
    ----------
    flood_results : pandas.DataFrame or geopandas.GeoDataFrame
        Flood hazard results returned by run_flood_hazard().

    output_dir : str or pathlib.Path
        Folder where flood-frequency plots will be saved.

    x_max : float, optional
        Fixed maximum x-axis value in expected events/year.
        If None, calculated automatically from all results.

    y_max : float, optional
        Fixed maximum y-axis value in expected events/year.
        If None, calculated automatically from all results.

    axis_padding : float, default 0.05
        Additional fraction added to automatically calculated axis limits.
        Example:
            0.05 = 5% padding.

    dpi : int, default 300
        Resolution of saved PNG files.

    show : bool, default False
        If True, display each plot.
        If False, save and close each plot.

    Returns
    -------
    list[pathlib.Path]
        Paths to created plot files.
    """

    # =========================================================================
    # IMPORTS
    # =========================================================================

    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd


    # =========================================================================
    # PREPARE OUTPUT DIRECTORY
    # =========================================================================

    output_dir = Path(
        output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    # =========================================================================
    # COPY RESULTS
    #
    # Work on a copy so plotting does not modify the original dataframe.
    # =========================================================================

    df = flood_results.copy()


    # =========================================================================
    # CHECK REQUIRED COLUMNS
    # =========================================================================

    required_columns = [
        "point_id",
        "name",
        "context",
        "duration_hr",
        "current_ari_year",
        "current_event_rate_per_year",
        "scenario",
        "time_horizon",
        "future_event_rate_per_year",
    ]


    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]


    if missing_columns:

        raise ValueError(
            "Flood plotting cannot continue because the following "
            "required columns are missing:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_columns
            )
        )


    # =========================================================================
    # AUTOMATICALLY DETECT CURRENT ARIs
    #
    # Do NOT hard-code 1, 10, 100.
    #
    # If the workflow uses:
    #
    #     [5, 10, 100]
    #
    # this automatically becomes:
    #
    #     [100, 10, 5]
    #
    # for plotting and labelling.
    # =========================================================================

    reference_aris = sorted(
        pd.to_numeric(
            df[
                "current_ari_year"
            ],
            errors="coerce",
        )
        .dropna()
        .unique(),
        reverse=True,
    )


    print(
        "[FLOOD PLOT] "
        "Reference ARIs detected: "
        f"{reference_aris}"
    )


    # =========================================================================
    # FUTURE RESULTS ONLY
    #
    # Historical rows are useful as references but should not generate
    # separate scenario/horizon figures.
    # =========================================================================

    future_df = df[
        df[
            "scenario"
        ]
        .astype(str)
        .str.lower()
        != "historical"
    ].copy()


    if future_df.empty:

        print(
            "[FLOOD PLOT] "
            "No future flood results available."
        )

        return []


    # =========================================================================
    # CONVERT FREQUENCY FIELDS TO NUMERIC
    # =========================================================================

    df[
        "current_event_rate_per_year"
    ] = pd.to_numeric(
        df[
            "current_event_rate_per_year"
        ],
        errors="coerce",
    )


    future_df[
        "future_event_rate_per_year"
    ] = pd.to_numeric(
        future_df[
            "future_event_rate_per_year"
        ],
        errors="coerce",
    )


    # =========================================================================
    # GLOBAL AXIS LIMITS
    #
    # IMPORTANT:
    #
    # Calculate x and y limits ONCE from the complete analysis.
    #
    # Do not calculate them separately for each plot.
    #
    # This means all SSP/horizon plots have identical scales and can be
    # compared visually.
    # =========================================================================

    current_values = df[
        "current_event_rate_per_year"
    ].to_numpy(
        dtype=float
    )


    future_values = future_df[
        "future_event_rate_per_year"
    ].to_numpy(
        dtype=float
    )


    # -------------------------------------------------------------------------
    # Remove NaN and infinite values.
    # -------------------------------------------------------------------------

    current_values = current_values[
        np.isfinite(
            current_values
        )
    ]


    future_values = future_values[
        np.isfinite(
            future_values
        )
    ]


    # -------------------------------------------------------------------------
    # AUTOMATIC X MAXIMUM
    #
    # x-axis = current expected number of events/year.
    # -------------------------------------------------------------------------

    if x_max is None:

        if current_values.size > 0:

            x_data_max = float(
                np.max(
                    current_values
                )
            )


            x_max_plot = (
                x_data_max
                * (
                    1.0
                    + axis_padding
                )
            )

        else:

            x_max_plot = 1.0


    else:

        x_max_plot = float(
            x_max
        )


    # -------------------------------------------------------------------------
    # AUTOMATIC Y MAXIMUM
    #
    # y-axis = future expected number of events/year.
    #
    # Include current values as well so the 1:1 reference line has a
    # sensible plotting extent.
    # -------------------------------------------------------------------------

    if y_max is None:

        all_frequency_values = np.concatenate(
            [
                current_values,
                future_values,
            ]
        )


        if all_frequency_values.size > 0:

            y_data_max = float(
                np.max(
                    all_frequency_values
                )
            )


            y_max_plot = (
                y_data_max
                * (
                    1.0
                    + axis_padding
                )
            )

        else:

            y_max_plot = 1.0


    else:

        y_max_plot = float(
            y_max
        )


    # =========================================================================
    # SAFETY CHECKS FOR AXES
    # =========================================================================

    if (
        not np.isfinite(
            x_max_plot
        )
        or x_max_plot <= 0
    ):

        x_max_plot = 1.0


    if (
        not np.isfinite(
            y_max_plot
        )
        or y_max_plot <= 0
    ):

        y_max_plot = 1.0


    print(
        "[FLOOD PLOT] "
        f"Common x-axis range: "
        f"0 to {x_max_plot:.4f} events/year"
    )


    print(
        "[FLOOD PLOT] "
        f"Common y-axis range: "
        f"0 to {y_max_plot:.4f} events/year"
    )


    # =========================================================================
    # GET UNIQUE SCENARIO / HORIZON COMBINATIONS
    # =========================================================================

    scenario_horizons = (
        future_df[
            [
                "scenario",
                "time_horizon",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "scenario",
                "time_horizon",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    created_plots = []


    # =========================================================================
    # LOOP THROUGH SCENARIOS / HORIZONS
    # =========================================================================

    for _, combination in scenario_horizons.iterrows():

        scenario = combination[
            "scenario"
        ]

        time_horizon = combination[
            "time_horizon"
        ]


        # ---------------------------------------------------------------------
        # Subset this scenario/horizon.
        # ---------------------------------------------------------------------

        plot_df = future_df[
            (
                future_df[
                    "scenario"
                ]
                == scenario
            )
            &
            (
                future_df[
                    "time_horizon"
                ]
                == time_horizon
            )
        ].copy()


        if plot_df.empty:

            continue


        # =====================================================================
        # CREATE FIGURE
        # =====================================================================

        fig, ax = plt.subplots(
            figsize=(
                9,
                6,
            )
        )


        # =====================================================================
        # GROUP BY INDIVIDUAL LOCATION / CONTEXT / DURATION
        #
        # This avoids averaging multiple locations.
        #
        # If later you add several urban or rural points, each one will be
        # represented separately.
        # =====================================================================

        group_columns = [
            "point_id",
            "name",
            "context",
            "duration_hr",
        ]


        for (
            point_id,
            point_name,
            context,
            duration_hr,
        ), group_df in plot_df.groupby(
            group_columns,
            dropna=False,
        ):


            # -----------------------------------------------------------------
            # Keep only the ARIs that were actually analysed.
            # -----------------------------------------------------------------

            group_df = group_df[
                group_df[
                    "current_ari_year"
                ].isin(
                    reference_aris
                )
            ].copy()


            # -----------------------------------------------------------------
            # Ensure numeric values.
            # -----------------------------------------------------------------

            group_df[
                "current_event_rate_per_year"
            ] = pd.to_numeric(
                group_df[
                    "current_event_rate_per_year"
                ],
                errors="coerce",
            )


            group_df[
                "future_event_rate_per_year"
            ] = pd.to_numeric(
                group_df[
                    "future_event_rate_per_year"
                ],
                errors="coerce",
            )


            # -----------------------------------------------------------------
            # Remove invalid rows.
            # -----------------------------------------------------------------

            valid = (
                np.isfinite(
                    group_df[
                        "current_event_rate_per_year"
                    ]
                )
                &
                np.isfinite(
                    group_df[
                        "future_event_rate_per_year"
                    ]
                )
            )


            group_df = group_df[
                valid
            ].copy()


            if group_df.empty:

                continue


            # -----------------------------------------------------------------
            # Sort by current event frequency.
            #
            # Example for:
            #
            #     100-year -> 0.01 EY
            #      10-year -> 0.11 EY
            #       5-year -> 0.20 EY
            #
            # left-to-right order becomes:
            #
            #     100-year
            #     10-year
            #     5-year
            # -----------------------------------------------------------------

            group_df = group_df.sort_values(
                "current_event_rate_per_year"
            )


            # =================================================================
            # CREATE LEGEND LABEL
            # =================================================================

            context_text = str(
                context
            ).strip().lower()


            if context_text == "urban":

                context_label = (
                    "urban losses"
                )


            elif context_text == "rural":

                context_label = (
                    "rural losses"
                )


            else:

                context_label = (
                    f"{context_text} losses"
                )


            # -----------------------------------------------------------------
            # If there are only two locations, the old-style label is fine:
            #
            #     1 hour event duration, urban losses
            #
            # If you later have multiple locations of the same context,
            # including the point name in the label is more useful.
            # -----------------------------------------------------------------

            same_context_count = (
                plot_df[
                    plot_df[
                        "context"
                    ]
                    == context
                ][
                    "point_id"
                ]
                .nunique()
            )


            if same_context_count > 1:

                label = (
                    f"{point_name}: "
                    f"{duration_hr:g} hour event duration, "
                    f"{context_label}"
                )


            else:

                label = (
                    f"{duration_hr:g} hour event duration, "
                    f"{context_label}"
                )


            # =================================================================
            # PLOT CURRENT VS FUTURE FREQUENCY
            # =================================================================

            ax.plot(
                group_df[
                    "current_event_rate_per_year"
                ],
                group_df[
                    "future_event_rate_per_year"
                ],
                marker="o",
                linewidth=2,
                markersize=7,
                label=label,
            )


            # =================================================================
            # LABEL EACH POINT WITH CURRENT ARI
            #
            # Examples:
            #
            #     5-yr
            #     10-yr
            #     100-yr
            # =================================================================

            for _, row in group_df.iterrows():

                ari = float(
                    row[
                        "current_ari_year"
                    ]
                )


                x_value = float(
                    row[
                        "current_event_rate_per_year"
                    ]
                )


                y_value = float(
                    row[
                        "future_event_rate_per_year"
                    ]
                )


                ax.annotate(
                    f"{ari:g}-yr",
                    xy=(
                        x_value,
                        y_value,
                    ),
                    xytext=(
                        6,
                        6,
                    ),
                    textcoords="offset points",
                    fontsize=8,
                )


        # =====================================================================
        # 1:1 NO-CHANGE LINE
        #
        # y = x
        #
        # Above this line:
        #     future event frequency is greater
        #
        # Below this line:
        #     future event frequency is smaller
        # =====================================================================

        one_to_one_max = min(
            x_max_plot,
            y_max_plot,
        )


        ax.plot(
            [
                0,
                one_to_one_max,
            ],
            [
                0,
                one_to_one_max,
            ],
            linestyle="--",
            linewidth=1.5,
            label="No change (1:1)",
        )


        # =====================================================================
        # COMMON AXIS LIMITS
        #
        # All figures use exactly these same limits.
        # =====================================================================

        ax.set_xlim(
            0,
            x_max_plot,
        )


        ax.set_ylim(
            0,
            y_max_plot,
        )


        # =====================================================================
        # AXIS LABELS
        # =====================================================================

        ax.set_xlabel(
            (
                "Expected number of events per year "
                "(current)"
            ),
            fontsize=12,
        )


        # ---------------------------------------------------------------------
        # Make scenario labels easier to read.
        #
        # Example:
        #
        #     ssp245 -> SSP2-4.5
        #     ssp370 -> SSP3-7.0
        # ---------------------------------------------------------------------

        scenario_labels = {
            "ssp126": "SSP1-2.6",
            "ssp245": "SSP2-4.5",
            "ssp370": "SSP3-7.0",
            "ssp585": "SSP5-8.5",
        }


        scenario_label = scenario_labels.get(
            str(
                scenario
            ).lower(),
            str(
                scenario
            ).upper(),
        )


        # ---------------------------------------------------------------------
        # Horizon label.
        # ---------------------------------------------------------------------

        horizon_label = str(
            time_horizon
        ).replace(
            "_",
            " ",
        ).title()

        ax.set_ylabel(
            (
                "Expected number of events per year "
                f"({scenario_label})"
            ),
            fontsize=12,
        )
        # =====================================================================
        # TITLE
        #
        # Read the current and future periods directly from TIME_WINDOWS.
        # This avoids hard-coding years in the plotting function.
        # =====================================================================

        baseline_start, baseline_end = time_windows[
            "baseline"
        ]

        future_start, future_end = time_windows[
            time_horizon
        ]


        ax.set_title(
            (
                f"Current ({baseline_start}–{baseline_end}) and future "
                f"({future_start}–{future_end}) flood frequencies"
            ),
            fontsize=12, fontweight="bold",
        )

        # =====================================================================
        # GRID
        # =====================================================================

        ax.grid(
            True,
            alpha=0.3,
        )


        # =====================================================================
        # LEGEND
        # =====================================================================

        handles, labels = ax.get_legend_handles_labels()


        if handles:

            ax.legend(
                loc="best",
                fontsize=10,
            )


        # =====================================================================
        # SAVE OUTPUT
        # =====================================================================

        safe_scenario = str(
            scenario
        ).replace(
            " ",
            "_",
        )


        safe_horizon = str(
            time_horizon
        ).replace(
            " ",
            "_",
        )


        output_path = (
            output_dir
            /
            (
                "Flood_Frequency_"
                f"{safe_scenario}_"
                f"{safe_horizon}.png"
            )
        )


        fig.tight_layout()


        fig.savefig(
            output_path,
            dpi=dpi,
            bbox_inches="tight",
        )


        # =====================================================================
        # SHOW OR CLOSE
        # =====================================================================

        if show:

            plt.show()


        else:

            plt.close(
                fig
            )


        # =====================================================================
        # RECORD OUTPUT
        # =====================================================================

        created_plots.append(
            output_path
        )


        print(
            "[FLOOD PLOT] "
            f"{output_path}"
        )


    # =========================================================================
    # FINISHED
    # =========================================================================

    print(
        f"Created {len(created_plots)} "
        "flood-frequency plot(s)."
    )


    return created_plots


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def run_flood_hazard(
    *,
    output_root: str | Path,
    time_windows: dict[str, tuple[int, int]],
    scenarios: list[str],

    # -----------------------------------------------------------------
    # LOCATION INPUT
    #
    # Two options are supported:
    #
    # 1. coordinate_points
    #       User-defined latitude/longitude locations.
    #
    # 2. location_path
    #       Point shapefile or GeoPackage.
    #
    # If coordinate_points is supplied, it takes priority.
    # Otherwise location_path is used.
    # -----------------------------------------------------------------
    coordinate_points: list[dict] | None = None,
    location_path: str | Path | None = None,

    name_field: str = "Name",
    context_field: str = "Name",
    urban_durations_hours: list[float] | None = None,
    rural_durations_hours: list[float] | None = None,
    current_aris_years: list[float] | None = None,
    timesteps: int = 60,
    urban_initial_loss_mm: float | None = None,
    urban_continuing_loss_mm_per_hr: float | None = None,
    timeout: int = 120,
    overwrite: bool = False,
    verbose: bool = True,
) -> gpd.GeoDataFrame:
    """
    Run the approximate climate-adjusted flood-frequency workflow.

    The method compares the runoff generated by a current design rainfall
    event with runoff generated by future climate-adjusted rainfall events
    across a broad range of frequencies.

    The future event frequency that produces approximately the same runoff
    as the current event is reported as the equivalent future flood
    frequency.

    ----------------------------------------------------------------------
    IMPORTANT
    ----------------------------------------------------------------------

    This is NOT a hydraulic flood model.

    It does not calculate:

        - flood extent;
        - flood depth;
        - flood velocity;
        - channel routing;
        - stormwater-network behaviour;
        - culvert capacity;
        - floodplain storage;
        - property-scale inundation.

    The output is an approximate climate-driven change in equivalent
    rainfall-runoff event frequency.

    ----------------------------------------------------------------------
    INPUT LOCATIONS
    ----------------------------------------------------------------------

    ``location_path`` must contain point features.

    Any number of points can be supplied.

    ``context_field`` must identify each point as Urban or Rural.

    For the current ACT example:

        Name
        -----
        Urban
        Rural

    ----------------------------------------------------------------------
    DEFAULT DURATIONS
    ----------------------------------------------------------------------

    Urban:
        1 hour

    Rural:
        24 hours

    These defaults reproduce the flood methodology's simplified national methodology but
    can be changed by the user.

    Example:

        urban_durations_hours=[
            1.0,
            4.5,
            12.0,
        ]

        rural_durations_hours=[
            24.0,
            72.0,
        ]

    ----------------------------------------------------------------------
    CURRENT REFERENCE EVENTS
    ----------------------------------------------------------------------

    Default ARIs:

        1 year
        10 years
        100 years

    These can also be changed by the user.

    ----------------------------------------------------------------------
    BOM DESIGN RAINFALL
    ----------------------------------------------------------------------

    Three BoM design-rainfall ranges are downloaded:

        very_frequent
        ifds
        rare

    These are combined to provide a broad rainfall-frequency curve for
    future runoff matching.

    ----------------------------------------------------------------------
    ARR INPUTS
    ----------------------------------------------------------------------

    ARR Data Hub provides:

        - current storm initial loss;
        - current storm continuing loss;
        - temperature change by horizon/scenario;
        - per-degree initial-loss adjustment factors;
        - per-degree continuing-loss adjustment factors.

    the flood methodology's urban calculation uses fixed losses of IL=2 mm and CL=0 mm/h.
    These are now the defaults. Optional urban override arguments are retained
    only so the workflow can be sensitivity-tested later.

    Rural losses are taken from ARR and climate-adjusted using the flood methodology's
    ``current_loss * (per_degree_factor ** deltaT)`` formulation.

    ----------------------------------------------------------------------
    OUTPUTS
    ----------------------------------------------------------------------

    Downloaded / source data:

        <output_root>/
            NARCliM Data/
                Flood/

    Final flood hazard outputs:

        <output_root>/
            Climate_Indicies/
                Hazards/
                    Flood/

    Outputs include:

        Flood_Hazard.gpkg
        Flood_Hazard.csv
        Flood_Runoff_Curves.csv

    Each point also receives QA/source files such as:

        ARR_Data_Hub.html
        ARR_losses_parsed.csv

        BoM_Very_Frequent.html
        BoM_IFD.html
        BoM_Rare.html
        BoM_Design_Rainfall_Combined.csv
    """

    # ==================================================================
    # DEFAULT USER SETTINGS
    # ==================================================================

    if urban_durations_hours is None:

        urban_durations_hours = [
            1.0,
        ]

    if rural_durations_hours is None:

        rural_durations_hours = [
            24.0,
        ]

    if current_aris_years is None:

        current_aris_years = [
            1.0,
            10.0,
            100.0,
        ]

    # ==================================================================
    # OUTPUT FOLDERS
    #
    # These folders are created even when the user runs ONLY flood and
    # has not previously run workflow.py or uncertainty.py.
    # ==================================================================

    output_root = Path(
        output_root
    )

    input_folder = (
        output_root
        / "NARCliM Data"
        / "Flood"
    )

    hazard_folder = (
        output_root
        / "Climate_Indicies"
        / "Hazards"
        / "Flood"
    )

    input_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    hazard_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ==================================================================
    # START LOGGING
    # ==================================================================

    _log(
        "=" * 100,
        verbose=verbose,
    )

    _log(
        "FLOOD HAZARD WORKFLOW",
        verbose=verbose,
    )

    _log(
        "=" * 100,
        verbose=verbose,
    )

    _log(
        f"[CODE VERSION] {FLOOD_CODE_VERSION}",
        verbose=verbose,
    )

    _log(
        f"[INPUT DATA] {input_folder}",
        verbose=verbose,
    )

    _log(
        f"[HAZARD OUTPUT] {hazard_folder}",
        verbose=verbose,
    )

    # ==================================================================
    # READ FLOOD LOCATIONS
    #
    # Priority:
    #     1. User-defined coordinates
    #     2. Shapefile / GeoPackage
    # ==================================================================

    locations, points = (
        _load_flood_points(
            location_path=location_path,
            coordinate_points=coordinate_points,
            name_field=name_field,
            context_field=context_field,
        )
    )

    if coordinate_points is not None:

        _log(
            "[LOCATION INPUT] "
            "Using user-defined latitude/longitude coordinates.",
            verbose=verbose,
        )

    else:

        _log(
            f"[LOCATION INPUT] "
            f"Using spatial file: {location_path}",
            verbose=verbose,
        )


    _log(
        f"[POINTS] {len(points)}",
        verbose=verbose,
    )

    _log(
        f"[POINTS] {len(points)}",
        verbose=verbose,
    )

    # ==================================================================
    # VALIDATE FUTURE SCENARIOS
    # ==================================================================

    for scenario in scenarios:

        if (
            scenario != "historical"
            and scenario
            not in ARR_SCENARIO_NAMES
        ):

            raise ValueError(
                "Unsupported ARR future scenario: "
                f"{scenario!r}. "
                "Supported future scenarios are "
                f"{sorted(ARR_SCENARIO_NAMES)}."
            )

    # ==================================================================
    # OUTPUT RECORD CONTAINERS
    # ==================================================================

    results: list[dict] = []

    runoff_curves: list[dict] = []

    # ==================================================================
    # PROCESS EACH FLOOD LOCATION
    # ==================================================================

    for point in points:

        _log(
            "\n"
            + "-" * 100,
            verbose=verbose,
        )

        _log(
            f"[POINT] "
            f"{point.name} "
            f"({point.context})",
            verbose=verbose,
        )

        # --------------------------------------------------------------
        # Create a source-data folder for this point.
        #
        # Example:
        #
        # NARCliM Data/
        #     Flood/
        #         F1_Urban/
        #         F2_Rural/
        # --------------------------------------------------------------

        point_input_folder = (
            input_folder
            / point.point_id
        )

        point_input_folder.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ==================================================================
        # ARR DATA HUB
        # ==================================================================

        arr_path = (
            point_input_folder
            / "ARR_Data_Hub.txt"
        )

        if (
            arr_path.exists()
            and not overwrite
        ):

            arr_html = (
                arr_path.read_text(
                    encoding="utf-8"
                )
            )

            arr_url = (
                "reused local "
                "ARR_Data_Hub.txt"
            )

            _log(
                "[ARR] Reusing saved ARR input.",
                verbose=verbose,
            )

        else:

            arr_html, arr_url = (
                _download_arr_html(
                    longitude=(
                        point.longitude
                    ),
                    latitude=(
                        point.latitude
                    ),
                    timeout=timeout,
                )
            )

            arr_path.write_text(
                arr_html,
                encoding="utf-8",
            )

            _log(
                "[ARR] Downloaded ARR input.",
                verbose=verbose,
            )

        # ==================================================================
        # BOM DESIGN RAINFALL
        #
        # the flood methodology's methodology requires a broad future-frequency search.
        #
        # Therefore we download:
        #
        #   1. Very Frequent
        #   2. Standard IFD
        #   3. Rare
        #
        # rather than using only the normal IFD table.
        # ==================================================================

        bom_files = {

            "very_frequent":
                (
                    point_input_folder
                    / "BoM_Very_Frequent.html"
                ),

            "ifds":
                (
                    point_input_folder
                    / "BoM_IFD.html"
                ),

            "rare":
                (
                    point_input_folder
                    / "BoM_Rare.html"
                ),
        }

        bom_html_by_design = {}

        # --------------------------------------------------------------
        # Download or reuse each BoM rainfall-frequency product.
        # --------------------------------------------------------------

        for (
            design_range,
            html_path,
        ) in bom_files.items():

            if (
                html_path.exists()
                and not overwrite
            ):

                html = (
                    html_path.read_text(
                        encoding="utf-8"
                    )
                )

                source_url = (
                    "reused local "
                    f"{html_path.name}"
                )

                _log(
                    f"[BoM] Reusing "
                    f"{design_range} input.",
                    verbose=verbose,
                )

            else:

                html, source_url = (
                    _download_bom_design_rainfall_html(
                        longitude=(
                            point.longitude
                        ),
                        latitude=(
                            point.latitude
                        ),
                        design=(
                            design_range
                        ),
                        timeout=timeout,
                    )
                )

                html_path.write_text(
                    html,
                    encoding="utf-8",
                )

                _log(
                    f"[BoM] Downloaded "
                    f"{design_range} input.",
                    verbose=verbose,
                )

            bom_html_by_design[
                design_range
            ] = (
                html,
                source_url,
            )

        # ==================================================================
        # PARSE + COMBINE BOM RAINFALL PRODUCTS
        #
        # The resulting table should contain a continuous set of
        # duration/frequency/rainfall-depth relationships across:
        #
        # Very Frequent + IFD + Rare.
        # ==================================================================

        bom_ifd = (
            _parse_all_bom_design_rainfall(
                bom_html_by_design
            )
        )

        # Save combined BoM source data for QA and reproducibility.

        bom_combined_csv = (
            point_input_folder
            / "BoM_Design_Rainfall_Combined.csv"
        )

        bom_ifd.to_csv(
            bom_combined_csv,
            index=False,
        )

        _log(
            "[BoM] Combined design-rainfall "
            f"records: {len(bom_ifd)}",
            verbose=verbose,
        )

        # ==================================================================
        # PARSE ARR DATA
        # ==================================================================

        arr_losses = (
            _parse_arr_storm_losses(
                arr_html
            )
        )

        # the flood methodology obtains temperature changes and per-degree loss factors from
        # the ARR text sections, then applies the equations explicitly.
        temperature_table = _parse_arr_temperature_changes(arr_html)

        (
            il_table,
            cl_table,
        ) = _parse_arr_loss_tables(arr_html)

        # ==================================================================
        # SAVE PARSED ARR LOSSES
        # ==================================================================

        pd.DataFrame(
            [
                {
                    "point_id":
                        point.point_id,

                    "name":
                        point.name,

                    "context":
                        point.context,

                    "arr_initial_loss_mm":
                        arr_losses.initial_loss_mm,

                    "arr_continuing_loss_mm_per_hr":
                        (
                            arr_losses
                            .continuing_loss_mm_per_hr
                        ),

                    "arr_source_url":
                        arr_url,

                    "bom_very_frequent_source":
                        (
                            bom_html_by_design[
                                "very_frequent"
                            ][1]
                        ),

                    "bom_ifd_source":
                        (
                            bom_html_by_design[
                                "ifds"
                            ][1]
                        ),

                    "bom_rare_source":
                        (
                            bom_html_by_design[
                                "rare"
                            ][1]
                        ),
                }
            ]
        ).to_csv(
            point_input_folder
            / "ARR_losses_parsed.csv",
            index=False,
        )

        # ==================================================================
        # SELECT CURRENT LOSSES BASED ON URBAN / RURAL CONTEXT
        # ==================================================================

        (
            current_losses,
            loss_status,
        ) = (
            _current_losses_for_context(
                context=(
                    point.context
                ),
                arr_losses=(
                    arr_losses
                ),
                urban_initial_loss_mm=(
                    urban_initial_loss_mm
                ),
                urban_continuing_loss_mm_per_hr=(
                    urban_continuing_loss_mm_per_hr
                ),
            )
        )

        # --------------------------------------------------------------
        # ARR storm losses are primarily rural.
        #
        # Allow the code to run for QA when urban overrides have not yet
        # been supplied, but clearly flag this situation.
        # --------------------------------------------------------------

        if (
            loss_status
            == "urban_using_arr_rural_losses_WARNING"
        ):

            warnings.warn(
                (
                    f"Urban point {point.name!r} "
                    "is currently using ARR rural losses "
                    "because no urban loss override was supplied. "
                    "Replace these with locally justified urban losses "
                    "before treating the urban flood result as final."
                ),
                RuntimeWarning,
            )

        # ==================================================================
        # SELECT DURATIONS
        # ==================================================================

        if (
            point.context
            == "urban"
        ):

            durations = (
                urban_durations_hours
            )

        else:

            durations = (
                rural_durations_hours
            )

        # ==================================================================
        # PROCESS EACH DURATION
        # ==================================================================

        for duration_hr in durations:

            duration_hr = float(
                duration_hr
            )

            # --------------------------------------------------------------
            # Extract the complete BoM frequency curve for this duration.
            #
            # Example:
            #
            # Urban:
            #     1 hour
            #
            # Rural:
            #     24 hours
            #
            # but the user can supply any duration available in BoM.
            # --------------------------------------------------------------

            duration_ifd = (
                _ifd_for_duration(
                    bom_ifd,
                    duration_hr,
                )
            )

            _log(
                f"[DURATION] {duration_hr:g} h | "
                f"{len(duration_ifd)} "
                "BoM frequency point(s)",
                verbose=verbose,
            )

            # ==================================================================
            # CURRENT REFERENCE EVENTS
            # ==================================================================

            for current_ari in (
                current_aris_years
            ):

                current_ari = float(
                    current_ari
                )

                if (
                    current_ari
                    <= 0
                ):

                    raise ValueError(
                        "Current ARIs must be > 0."
                    )
                # ----------------------------------------------------------
                # reference event frequency.
                #
                # 1-year   -> 1.00 EY
                # 10-year  -> 0.11 EY
                # 100-year -> 0.01 EY
                # ----------------------------------------------------------
                current_rate = _reference_event_rate(current_ari)

                # ----------------------------------------------------------
                # Interpolate current BoM rainfall depth for this duration
                # and reference frequency.
                # ----------------------------------------------------------

                current_depth = (
                    _interpolate_ifd_depth(
                        duration_ifd,
                        current_rate,
                    )
                )

                # ----------------------------------------------------------
                # Calculate CURRENT rainfall excess/runoff.
                #
                # Rainfall is distributed across the user-defined number
                # of timesteps (default 60).
                #
                # Initial loss is removed first, followed by continuing
                # loss. Remaining rainfall becomes runoff.
                # ----------------------------------------------------------

                current_runoff = (
                    _total_runoff(
                        rainfall_depth_mm=(
                            current_depth
                        ),
                        duration_hours=(
                            duration_hr
                        ),
                        initial_loss_mm=(
                            current_losses
                            .initial_loss_mm
                        ),
                        continuing_loss_mm_per_hr=(
                            current_losses
                            .continuing_loss_mm_per_hr
                        ),
                        timesteps=timesteps,
                    )
                )
                # ==============================================================
                # CURRENT RUNOFF THRESHOLD / RAINFALL SHORTFALL
                # rainfall excess is >= 0.1 mm.  For lower-runoff events it
                # attempts a rainfall-shortfall fallback.  Here we calculate
                # that shortfall robustly as the extra total storm depth needed
                # to produce 0.1 mm runoff using the SAME hydrograph routine.
                # ==============================================================

                (
                    current_rainfall_shortfall_mm,
                    current_runoff_threshold_depth_mm,
                ) = _rainfall_shortfall_to_runoff(
                    rainfall_depth_mm=current_depth,
                    duration_hours=duration_hr,
                    initial_loss_mm=current_losses.initial_loss_mm,
                    continuing_loss_mm_per_hr=(
                        current_losses.continuing_loss_mm_per_hr
                    ),
                    timesteps=timesteps,
                )

                current_runoff_is_zero = (
                    not np.isfinite(current_runoff)
                    or current_runoff < RUNOFF_MATCH_THRESHOLD_MM
                )
                # ==================================================================
                # TIME HORIZONS
                # ==================================================================

                for (
                    horizon,
                    (
                        start_year,
                        end_year,
                    ),
                ) in (
                    time_windows.items()
                ):

                    # ----------------------------------------------------------
                    # Use the representative/midpoint year to obtain ARR
                    # climate-change factors.
                    #
                    # Examples:
                    #
                    # 1985-2014 -> 2000
                    # 2041-2060 -> 2050
                    # 2081-2100 -> 2090
                    # ----------------------------------------------------------

                    mid_year = (
                        _representative_year(
                            start_year,
                            end_year,
                        )
                    )

                    # ==========================================================
                    # SCENARIOS
                    # ==========================================================

                    for scenario in scenarios:

                        # ------------------------------------------------------
                        # Historical scenario is paired only with historical
                        # time periods.
                        # ------------------------------------------------------

                        if (
                            scenario
                            == "historical"
                            and end_year
                            > 2014
                        ):

                            continue

                        # ------------------------------------------------------
                        # SSP scenarios are paired only with future periods.
                        # ------------------------------------------------------

                        if (
                            scenario
                            != "historical"
                            and end_year
                            <= 2014
                        ):

                            continue

                        # ------------------------------------------------------
                        # Fields shared by historical and future outputs.
                        # ------------------------------------------------------

                        common = {

                            "point_id":
                                point.point_id,

                            "name":
                                point.name,

                            "context":
                                point.context,

                            "duration_hr":
                                duration_hr,

                            "current_ari_year":
                                current_ari,

                            "current_event_rate_per_year":
                                current_rate,

                            "current_rainfall_depth_mm":
                                current_depth,

                            "current_initial_loss_mm":
                                (
                                    current_losses
                                    .initial_loss_mm
                                ),

                            "current_continuing_loss_mm_per_hr":
                                (
                                    current_losses
                                    .continuing_loss_mm_per_hr
                                ),

                            "current_runoff_mm":
                                current_runoff,

                            "current_runoff_threshold_mm":
                                RUNOFF_MATCH_THRESHOLD_MM,

                            "current_runoff_threshold_depth_mm":
                                current_runoff_threshold_depth_mm,

                            "current_rainfall_shortfall_mm":
                                current_rainfall_shortfall_mm,

                            "scenario":
                                scenario,

                            "time_horizon":
                                horizon,

                            "start_year":
                                start_year,

                            "end_year":
                                end_year,

                            "time_horizon_year":
                                mid_year,

                            "loss_source":
                                current_losses.source,

                            "loss_status":
                                loss_status,
                        }

                        # ======================================================
                        # HISTORICAL BASELINE
                        #
                        # No climate adjustment is applied.
                        #
                        # The equivalent future/current frequency is identical
                        # by definition.
                        # ======================================================

                        if (
                            scenario
                            == "historical"
                        ):

                            results.append(
                                {
                                    **common,

                                    "delta_t_degC":
                                        0.0,

                                    "rainfall_alpha_pct_per_degC":
                                        _rainfall_alpha_for_duration(duration_hr),

                                    "rainfall_factor":
                                        1.0,

                                    "initial_loss_per_degree_factor":
                                        1.0,

                                    "continuing_loss_per_degree_factor":
                                        1.0,

                                    "initial_loss_factor":
                                        1.0,

                                    "continuing_loss_factor":
                                        1.0,

                                    "future_initial_loss_mm":
                                        (
                                            current_losses
                                            .initial_loss_mm
                                        ),

                                    "future_continuing_loss_mm_per_hr":
                                        (
                                            current_losses
                                            .continuing_loss_mm_per_hr
                                        ),

                                    "future_runoff_threshold_depth_mm":
                                        current_runoff_threshold_depth_mm,

                                    "match_basis":
                                        "historical_reference",

                                    "future_event_rate_per_year":
                                        current_rate,

                                    "future_ari_year":
                                        current_ari,

                                    "frequency_ratio":
                                        1.0,

                                    "frequency_change_pct":
                                        0.0,

                                    "match_status":
                                        "historical_reference",
                                }
                            )

                            continue

                        # ======================================================
                        # FLOOD CLIMATE-CHANGE CALCULATIONS
                        # ======================================================

                        # Temperature change is read from ARR for this horizon and scenario.
                        delta_t = _arr_delta_t(
                            temperature_table,
                            year=mid_year,
                            scenario=scenario,
                        )

                        # the flood methodology rainfall adjustment:
                        # future rain = current rain * (1 + alpha/100) ** deltaT
                        rainfall_alpha_pct_per_deg = _rainfall_alpha_for_duration(
                            duration_hr
                        )
                        rainfall_factor = _rainfall_climate_factor(
                            duration_hours=duration_hr,
                            delta_t=delta_t,
                        )

                        if point.context == "urban":
                            # the flood methodology keeps urban losses fixed in the future:
                            # IL = 2 mm, CL = 0 mm/h (or supplied urban override).
                            il_per_degree_factor = 1.0
                            cl_per_degree_factor = 1.0
                            il_factor = 1.0
                            cl_factor = 1.0
                            future_il = current_losses.initial_loss_mm
                            future_cl = current_losses.continuing_loss_mm_per_hr
                        else:
                            # the flood methodology reads ARR loss adjustment factors and then
                            # compounds them by the temperature change:
                            # future loss = current loss * (per_degree_factor ** dT)
                            il_per_degree_factor = _arr_loss_per_degree_factor(
                                il_table,
                                year=mid_year,
                                scenario=scenario,
                            )
                            cl_per_degree_factor = _arr_loss_per_degree_factor(
                                cl_table,
                                year=mid_year,
                                scenario=scenario,
                            )

                            il_factor = il_per_degree_factor ** delta_t
                            cl_factor = cl_per_degree_factor ** delta_t

                            future_il = (
                                current_losses.initial_loss_mm
                                * il_factor
                            )
                            future_cl = (
                                current_losses.continuing_loss_mm_per_hr
                                * cl_factor
                            )

                        # ------------------------------------------------------
                        # Rainfall depth required to produce 0.1 mm future
                        # runoff under this horizon/scenario loss regime.
                        # This threshold is constant across all future EYs for
                        # the same point, duration, horizon and scenario.
                        # ------------------------------------------------------
                        future_runoff_threshold_depth_mm = (
                            _minimum_rainfall_depth_for_runoff(
                                duration_hours=duration_hr,
                                initial_loss_mm=future_il,
                                continuing_loss_mm_per_hr=future_cl,
                                timesteps=timesteps,
                            )
                        )

                        # ======================================================
                        # NEAR-ZERO CURRENT RUNOFF
                        #
                        # Do NOT return undefined here.  The revised workflow
                        # rainfall shortfall instead of runoff.
                        # ======================================================

                        if current_runoff_is_zero:
                            _log(
                                (
                                    "[QA] Current runoff is < 0.1 mm; "
                                    "using corrected rainfall-shortfall "
                                    "frequency matching. "
                                    f"Current shortfall="
                                    f"{current_rainfall_shortfall_mm:.3f} mm."
                                ),
                                verbose=verbose,
                            )


                        # ======================================================
                        # BUILD FUTURE RUNOFF CURVE
                        #
                        # For every available BoM frequency:
                        #
                        # current design rainfall depth
                        #       ↓
                        # climate rainfall factor
                        #       ↓
                        # future rainfall depth
                        #       ↓
                        # future IL / CL
                        #       ↓
                        # future runoff
                        #
                        # This produces the curve used to find the future
                        # equivalent event frequency.
                        # ======================================================

                        curve_rows = []

                        for (
                            _,
                            ifd_row,
                        ) in (
                            duration_ifd.iterrows()
                        ):

                            rate = float(
                                ifd_row[
                                    "event_rate_per_year"
                                ]
                            )

                            base_depth = float(
                                ifd_row[
                                    "rainfall_depth_mm"
                                ]
                            )

                            # --------------------------------------------------
                            # Climate-adjust design rainfall.
                            # --------------------------------------------------

                            future_depth = (
                                base_depth
                                * rainfall_factor
                            )

                            # --------------------------------------------------
                            # Calculate future rainfall excess/runoff.
                            # --------------------------------------------------

                            future_runoff = (
                                _total_runoff(
                                    rainfall_depth_mm=(
                                        future_depth
                                    ),
                                    duration_hours=(
                                        duration_hr
                                    ),
                                    initial_loss_mm=(
                                        future_il
                                    ),
                                    continuing_loss_mm_per_hr=(
                                        future_cl
                                    ),
                                    timesteps=(
                                        timesteps
                                    ),
                                )
                            )

                            # Additional total rainfall depth needed for this
                            # future event to reach the 0.1 mm runoff threshold.
                            future_rainfall_shortfall_mm = max(
                                0.0,
                                future_runoff_threshold_depth_mm
                                - future_depth,
                            )

                            row = {

                                "point_id":
                                    point.point_id,

                                "name":
                                    point.name,

                                "context":
                                    point.context,

                                "duration_hr":
                                    duration_hr,

                                "current_ari_year":
                                    current_ari,

                                "scenario":
                                    scenario,

                                "time_horizon":
                                    horizon,

                                "time_horizon_year":
                                    mid_year,

                                "event_rate_per_year":
                                    rate,

                                "current_ifd_depth_mm":
                                    base_depth,

                                "delta_t_degC":
                                    delta_t,

                                "rainfall_alpha_pct_per_degC":
                                    rainfall_alpha_pct_per_deg,

                                "rainfall_factor":
                                    rainfall_factor,

                                "initial_loss_per_degree_factor":
                                    il_per_degree_factor,

                                "continuing_loss_per_degree_factor":
                                    cl_per_degree_factor,

                                "future_rainfall_depth_mm":
                                    future_depth,

                                "future_initial_loss_mm":
                                    future_il,

                                "future_continuing_loss_mm_per_hr":
                                    future_cl,

                                "future_runoff_mm":
                                    future_runoff,

                                "future_runoff_threshold_depth_mm":
                                    future_runoff_threshold_depth_mm,

                                "future_rainfall_shortfall_mm":
                                    future_rainfall_shortfall_mm,
                            }

                            # Keep a local copy for interpolation.
                            curve_rows.append(
                                row
                            )

                            # Also retain the complete curve as a QA output.
                            runoff_curves.append(
                                row
                            )

                        future_curve = (
                            pd.DataFrame(
                                curve_rows
                            )
                        )

                        # ======================================================
                        # MATCH FUTURE FREQUENCY
                        #
                        # Solve approximately:
                        #
                        # future runoff(frequency)
                        #
                        #     =
                        #
                        # current runoff
                        #
                        # Interpolation is performed between available BoM
                        # event frequencies.
                        # ======================================================

                        (
                            future_rate,
                            match_status,
                        ) = (
                            _match_future_event_rate(
                                current_runoff_mm=(
                                    current_runoff
                                ),
                                current_rainfall_shortfall_mm=(
                                    current_rainfall_shortfall_mm
                                ),
                                future_curve=(
                                    future_curve
                                ),
                            )
                        )

                        # ======================================================
                        # CONVERT FUTURE EVENT RATE TO ARI
                        # ======================================================

                        if (
                            np.isfinite(
                                future_rate
                            )
                            and future_rate
                            > 0
                        ):

                            future_ari = (
                                1.0
                                / future_rate
                            )

                            # --------------------------------------------------
                            # Ratio > 1:
                            #
                            # equivalent event becomes MORE frequent.
                            #
                            # Example:
                            #
                            # current 1-in-10
                            # future 1-in-5
                            #
                            # current rate = 0.1
                            # future rate  = 0.2
                            #
                            # ratio = 2
                            #
                            # => 100% more frequent.
                            # --------------------------------------------------

                            ratio = (
                                future_rate
                                / current_rate
                            )

                            change_pct = (
                                (
                                    ratio
                                    - 1.0
                                )
                                * 100.0
                            )

                        else:

                            future_ari = (
                                np.nan
                            )

                            ratio = (
                                np.nan
                            )

                            change_pct = (
                                np.nan
                            )

                        # ======================================================
                        # SAVE FINAL FUTURE RECORD
                        # ======================================================

                        results.append(
                            {
                                **common,

                                "delta_t_degC":
                                    delta_t,

                                "rainfall_alpha_pct_per_degC":
                                    rainfall_alpha_pct_per_deg,

                                "rainfall_factor":
                                    rainfall_factor,

                                "initial_loss_per_degree_factor":
                                    il_per_degree_factor,

                                "continuing_loss_per_degree_factor":
                                    cl_per_degree_factor,

                                "initial_loss_factor":
                                    il_factor,

                                "continuing_loss_factor":
                                    cl_factor,

                                "future_initial_loss_mm":
                                    future_il,

                                "future_continuing_loss_mm_per_hr":
                                    future_cl,

                                "future_runoff_threshold_depth_mm":
                                    future_runoff_threshold_depth_mm,

                                "match_basis":
                                    (
                                        "rainfall_shortfall"
                                        if current_runoff_is_zero
                                        else "runoff"
                                    ),

                                "future_event_rate_per_year":
                                    future_rate,

                                "future_ari_year":
                                    future_ari,

                                "frequency_ratio":
                                    ratio,

                                "frequency_change_pct":
                                    change_pct,

                                "match_status":
                                    match_status,
                            }
                        )
    # ==================================================================
    # VALIDATE FINAL RESULTS
    # ==================================================================

    if not results:

        raise RuntimeError(
            "Flood workflow produced no results."
        )

    # ==================================================================
    # CREATE FINAL DATAFRAME
    # ==================================================================

    df = pd.DataFrame(
        results
    )

    # ==================================================================
    # ROUND FINAL OUTPUTS
    #
    # Keep enough precision for frequency calculations while keeping GIS
    # and CSV tables readable.
    # ==================================================================

    round_map = {

        "duration_hr":
            2,

        "current_ari_year":
            3,

        "current_event_rate_per_year":
            6,

        "current_rainfall_depth_mm":
            2,

        "current_initial_loss_mm":
            2,

        "current_continuing_loss_mm_per_hr":
            3,

        "current_runoff_mm":
            2,

        "delta_t_degC":
            3,

        "rainfall_alpha_pct_per_degC":
            3,

        "rainfall_factor":
            6,

        "initial_loss_per_degree_factor":
            6,

        "continuing_loss_per_degree_factor":
            6,

        "initial_loss_factor":
            3,

        "continuing_loss_factor":
            3,

        "future_initial_loss_mm":
            2,

        "future_continuing_loss_mm_per_hr":
            3,

        "future_event_rate_per_year":
            6,

        "future_ari_year":
            3,

        "frequency_ratio":
            3,

        "frequency_change_pct":
            1,
    }

    for (
        column,
        decimals,
    ) in (
        round_map.items()
    ):

        if (
            column
            in df.columns
        ):

            df[
                column
            ] = (
                pd.to_numeric(
                    df[
                        column
                    ],
                    errors="coerce",
                )
                .round(
                    decimals
                )
            )

    # ==================================================================
    # ADD POINT GEOMETRY
    # ==================================================================

    geometry_lookup = (
        locations[
            [
                "_flood_point_id",
                "geometry",
            ]
        ]
        .rename(
            columns={
                "_flood_point_id":
                    "point_id"
            }
        )
    )

    final = (
        gpd.GeoDataFrame(
            df.merge(
                geometry_lookup,
                on="point_id",
                how="left",
            ),
            geometry="geometry",
            crs="EPSG:4326",
        )
    )

    # ==================================================================
    # FINAL OUTPUT PATHS
    # ==================================================================

    output_gpkg = (
        hazard_folder
        / "Flood_Hazard.gpkg"
    )

    output_csv = (
        hazard_folder
        / "Flood_Hazard.csv"
    )

    runoff_curve_csv = (
        hazard_folder
        / "Flood_Runoff_Curves.csv"
    )

    # ==================================================================
    # WRITE GEOPACKAGE
    # ==================================================================

    if (
        output_gpkg.exists()
        and overwrite
    ):

        output_gpkg.unlink()

    final.to_file(
        output_gpkg,
        layer="flood_frequency_change",
        driver="GPKG",
        index=False,
    )

    # ==================================================================
    # WRITE CSV
    # ==================================================================

    pd.DataFrame(
        final.drop(
            columns="geometry"
        )
    ).to_csv(
        output_csv,
        index=False,
    )

    # ==================================================================
    # WRITE FULL RUNOFF CURVES FOR QA
    # ==================================================================

    runoff_curve_df = (
        pd.DataFrame(
            runoff_curves
        )
    )

    if (
        not runoff_curve_df.empty
    ):

        runoff_curve_df.to_csv(
            runoff_curve_csv,
            index=False,
        )

    # ==================================================================
    # FINISH
    # ==================================================================

    _log(
        "\n"
        + "=" * 100,
        verbose=verbose,
    )

    _log(
        "FLOOD HAZARD WORKFLOW FINISHED",
        verbose=verbose,
    )

    _log(
        "=" * 100,
        verbose=verbose,
    )

    _log(
        f"[GPKG] {output_gpkg}",
        verbose=verbose,
    )

    _log(
        f"[CSV] {output_csv}",
        verbose=verbose,
    )

    if (
        not runoff_curve_df.empty
    ):

        _log(
            f"[RUNOFF CURVES] "
            f"{runoff_curve_csv}",
            verbose=verbose,
        )

    return final