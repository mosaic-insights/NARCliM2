from __future__ import annotations

from pathlib import Path
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from matplotlib.patches import Patch


# =============================================================================
# DEFAULT CONFIGURATION
# =============================================================================

DEFAULT_CMAP = "viridis"
DEFAULT_DPI = 300

# Variable-specific units and labels. SPI12 is a dimensionless index.
VARIABLE_LABELS = {
    "r20mm": ("Days with rainfall >=20 mm", "days/year"),
    "txge35": ("Days with maximum temperature >=35°C", "days/year"),
    "ffdigt50": ("Days with FFDI >50", "days/year"),
    "stormdays": ("Storm days", "days/year"),
    "stormdaysgt89": ("Days with wind speed >89 km/h", "days/year"),
    "spi12": ("12-month Standardised Precipitation Index", "dimensionless"),
}

def _variable_info(raster_path: Path) -> tuple[str, str, str]:
    """Return variable, descriptive label and unit from the raster filename."""
    variable = re.split(r"_(?:Change_)?(?:Min|Max|Mean)_", raster_path.stem, maxsplit=1, flags=re.IGNORECASE)[0]
    label, unit = VARIABLE_LABELS.get(variable.lower(), (variable, "index units"))
    return variable, label, unit

def _is_change_raster(raster_path: Path) -> bool:
    return raster_path.parent.name.lower() == "change" or bool(
        re.search(r"_Change_(Min|Mean|Max)_", raster_path.stem, re.IGNORECASE)
    )

def _period_label(raster_path: Path) -> str:
    match = re.search(r"_(\d{4})_(\d{4})(?:_ID\d+)?$", raster_path.stem)
    if match:
        start = int(match.group(1))
        return "2050s" if 2035 <= start < 2060 else "2090s" if start >= 2070 else f"{match.group(1)}–{match.group(2)}"
    return "future projection"

def _overlay_layer(ax, path, raster_crs, *, color, label, alpha, linewidth):
    """Draw semi-transparent polygon fills and outlines; return a matching legend patch."""
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("overlay_alpha must be between 0 and 1")
    layer = gpd.read_file(path)
    if layer.empty:
        raise ValueError(f"Overlay contains no features: {path}")
    if layer.crs is None:
        raise ValueError(f"Overlay has no CRS: {path}")
    layer = layer.to_crs(raster_crs)
    geometry_types = set(layer.geometry.geom_type.dropna())
    if geometry_types <= {"Polygon", "MultiPolygon"}:
        # The alpha controls fill opacity; boundaries remain visible independently.
        layer.plot(ax=ax, facecolor=color, edgecolor="none", alpha=alpha, zorder=4)
        layer.boundary.plot(ax=ax, color=color, linewidth=linewidth, alpha=1.0, zorder=5)
        return Patch(facecolor=color, edgecolor=color, alpha=max(alpha, 0.35), label=label)
    if geometry_types <= {"LineString", "MultiLineString"}:
        layer.plot(ax=ax, color=color, linewidth=linewidth, zorder=5)
        from matplotlib.lines import Line2D
        return Line2D([0], [0], color=color, linewidth=linewidth, label=label)
    if geometry_types <= {"Point", "MultiPoint"}:
        layer.plot(ax=ax, color=color, markersize=max(4, linewidth * 12), zorder=5)
        from matplotlib.lines import Line2D
        return Line2D([0], [0], color=color, marker="o", linestyle="None", label=label)
    raise ValueError(f"Unsupported or mixed overlay geometry types {geometry_types}: {path}")


def _prepare_overlays(
    overlays: list[dict] | None,
    *,
    flood_extent_path: str | Path | None,
    bushfire_prone_path: str | Path | None,
    flood_color: str,
    bushfire_color: str,
    overlay_alpha: float,
) -> list[dict]:
    """Normalise user-defined overlays and preserve older flood/bushfire arguments.

    Each overlay accepts path (required), label (optional; defaults to filename),
    color, alpha, and linewidth. Each can be independently styled.
    """
    result = []
    for item in overlays or []:
        if not isinstance(item, dict) or not item.get("path"):
            raise ValueError("Each overlay must be a dictionary with a 'path'.")
        layer = dict(item)
        layer["path"] = Path(layer["path"])
        layer.setdefault("label", layer["path"].stem.replace("_", " "))
        layer.setdefault("color", "yellow")
        layer.setdefault("alpha", overlay_alpha)
        layer.setdefault("linewidth", 0.9)
        result.append(layer)

    # Backward compatibility with the previous two dedicated overlay options.
    if flood_extent_path is not None:
        result.append(dict(path=Path(flood_extent_path), label="Flood extent",
                           color=flood_color, alpha=overlay_alpha, linewidth=0.9))
    if bushfire_prone_path is not None:
        result.append(dict(path=Path(bushfire_prone_path), label="Bushfire-prone area",
                           color=bushfire_color, alpha=overlay_alpha, linewidth=0.9))
    return result


def _create_missing_spi12_change_rasters(ensemble_folder: Path, *, overwrite: bool = False, verbose: bool = True) -> None:
    """Create SPI12 projection-minus-historical GeoTIFFs if not already available.

    Match the same ensemble statistic (Min, Max, Mean) and feature ID.
    Reproject baseline to the future grid when grids differ.
    """
    rasters = sorted(ensemble_folder.rglob("SPI12_*.tif"))
    pattern = re.compile(r"^SPI12_(Min|Max|Mean)_(historical|ssp\d+)_(\d{4})_(\d{4})(_ID\d+)?$", re.IGNORECASE)
    baseline = {}
    futures = []
    for path in rasters:
        match = pattern.match(path.stem)
        if not match:
            continue
        stat, scenario, start, end, feature = match.groups()
        key = (stat.lower(), (feature or "").lower())
        if scenario.lower() == "historical":
            baseline[key] = path
        else:
            futures.append((path, stat, scenario, start, end, feature or "", key))
    if not futures:
        return
    change_folder = ensemble_folder / "Change"
    for future_path, stat, scenario, start, end, feature, key in futures:
        baseline_path = baseline.get(key)
        if baseline_path is None:
            if verbose:
                print(f"[SPI12 CHANGE SKIP] No historical {stat} baseline for {future_path.name}")
            continue
        output_path = change_folder / f"SPI12_Change_{stat}_{scenario}_{start}_{end}{feature}.tif"
        if output_path.exists() and not overwrite:
            continue
        with rasterio.open(future_path) as future, rasterio.open(baseline_path) as historical:
            if future.crs is None or historical.crs is None:
                raise ValueError("SPI12 baseline and projection rasters require CRS metadata")
            projected = future.read(1, masked=True).astype("float32")
            if (future.crs == historical.crs and future.transform == historical.transform
                    and future.shape == historical.shape):
                reference = historical.read(1, masked=True).astype("float32")
            else:
                aligned = np.full(future.shape, np.nan, dtype="float32")
                reproject(
                    source=rasterio.band(historical, 1), destination=aligned,
                    src_transform=historical.transform, src_crs=historical.crs,
                    src_nodata=historical.nodata, dst_transform=future.transform,
                    dst_crs=future.crs, dst_nodata=np.nan,
                    resampling=Resampling.bilinear,
                )
                reference = np.ma.masked_invalid(aligned)
            delta = np.ma.array(projected - reference)
            values = delta.filled(np.nan).astype("float32")
            values[~np.isfinite(values)] = np.nan
            profile = future.profile.copy()
            profile.update(dtype="float32", count=1, nodata=np.nan, compress="deflate")
            change_folder.mkdir(parents=True, exist_ok=True)
            with rasterio.open(output_path, "w", **profile) as dst:
                dst.write(values, 1)
                dst.update_tags(change_type="scenario_minus_baseline", baseline_raster=str(baseline_path),
                                projection_raster=str(future_path), variable="SPI12", statistic=stat)
        if verbose:
            print(f"[SPI12 CHANGE CREATED] {output_path.name}")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _safe_name(value: str) -> str:
    """Convert text into a safe filename component."""

    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value),
    ).strip("_")


def _raster_extent(
    src: rasterio.io.DatasetReader,
) -> tuple[float, float, float, float]:
    """Return raster extent as left, right, bottom, top for imshow()."""

    bounds = src.bounds
    return bounds.left, bounds.right, bounds.bottom, bounds.top


def _clean_raster_name(raster_path: Path) -> str:
    """
    Convert an Ensemble_Stats filename into a readable map title.

    Supported examples
    ------------------
    New naming convention:
        TXge35_Mean_ssp245_2040_2060.tif
        -> TXge35 | Mean | SSP245 | 2040-2060

    Change naming convention:
        TXge35_Change_Min_ssp245_2040_2059.tif
        -> TXge35 | Change in Min | SSP245 | 2040-2059

    Old naming convention:
        TXge35_ssp245_mid_term_ID1_ensemble_mean.tif
        -> TXge35 | SSP245 | Mid Term | Ensemble Mean
    """

    stem = raster_path.stem

    # ==================================================================
    # NEW NAMING CONVENTION
    # ==================================================================

    change_match = re.match(
        r"^(?P<variable>.+?)_Change_(?P<stat>Min|Mean|Max)_"
        r"(?P<scenario>historical|ssp[0-9]+)_"
        r"(?P<start_year>[0-9]{4})_(?P<end_year>[0-9]{4})"
        r"(?:_(?P<feature_id>.+))?$",
        stem,
        flags=re.IGNORECASE,
    )

    if change_match is not None:
        scenario = change_match.group("scenario")
        scenario = scenario.upper() if scenario.lower().startswith("ssp") else scenario.capitalize()
        return (
            f"{change_match.group('variable')} | "
            f"Change in {change_match.group('stat').capitalize()} | "
            f"{scenario} | "
            f"{change_match.group('start_year')}-{change_match.group('end_year')}"
        )

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

        variable = new_match.group(
            "variable"
        )

        statistic = new_match.group(
            "stat"
        ).capitalize()

        scenario = new_match.group(
            "scenario"
        )

        if scenario.lower().startswith(
            "ssp"
        ):
            scenario = scenario.upper()
        else:
            scenario = scenario.capitalize()

        start_year = new_match.group(
            "start_year"
        )

        end_year = new_match.group(
            "end_year"
        )

        return (
            f"{variable} | "
            f"{statistic} | "
            f"{scenario} | "
            f"{start_year}-{end_year}"
        )

    # ==================================================================
    # OLD NAMING CONVENTION
    # ==================================================================

    stem = re.sub(
        r"_ID\d+_",
        "_",
        stem,
        flags=re.IGNORECASE,
    )

    parts = stem.split("_")
    readable_parts: list[str] = []

    skip_next_term = False

    for index, part in enumerate(parts):

        if skip_next_term:
            skip_next_term = False
            continue

        lower = part.lower()

        if lower.startswith("ssp"):
            readable_parts.append(
                part.upper()
            )

        elif lower == "historical":
            readable_parts.append(
                "Historical"
            )

        elif lower == "baseline":
            readable_parts.append(
                "Baseline"
            )

        elif (
            lower in {"short", "mid", "long"}
            and index + 1 < len(parts)
            and parts[index + 1].lower() == "term"
        ):
            readable_parts.append(
                f"{part.capitalize()} Term"
            )
            skip_next_term = True

        elif lower == "ensemble":
            # Combine "ensemble" with the following statistic if present.
            if (
                index + 1 < len(parts)
                and parts[index + 1].lower()
                in {"mean", "min", "max"}
            ):
                stat = parts[
                    index + 1
                ].lower()

                stat_label = {
                    "mean": "Mean",
                    "min": "Minimum",
                    "max": "Maximum",
                }[
                    stat
                ]

                readable_parts.append(
                    f"Ensemble {stat_label}"
                )

                skip_next_term = True
            else:
                readable_parts.append(
                    "Ensemble"
                )

        elif lower == "mean":
            readable_parts.append(
                "Mean"
            )

        elif lower == "min":
            readable_parts.append(
                "Minimum"
            )

        elif lower == "max":
            readable_parts.append(
                "Maximum"
            )

        elif lower != "term":
            readable_parts.append(
                part
            )

    return " | ".join(
        readable_parts
    )

def _get_valid_values(raster: np.ma.MaskedArray) -> np.ndarray:
    """Return finite valid raster values."""

    values = raster.compressed()
    return values[np.isfinite(values)]


# =============================================================================
# PLOT ONE ENSEMBLE RASTER
# =============================================================================

def plot_ensemble_raster(
    *,
    raster_path: str | Path,
    boundary_path: str | Path,
    output_path: str | Path,
    cmap: str = DEFAULT_CMAP,
    dpi: int = DEFAULT_DPI,
    title: str | None = None,
    boundary_linewidth: float = 1.5,
    show: bool = False,
    flood_extent_path: str | Path | None = None,
    bushfire_prone_path: str | Path | None = None,
    flood_color: str = "lightskyblue",
    bushfire_color: str = "yellow",
    overlay_alpha: float = 0.22,
    overlays: list[dict] | None = None,
) -> Path:
    """
    Plot one NARCliM ensemble-statistic raster and overlay the study boundary.
    """

    raster_path = Path(raster_path)
    boundary_path = Path(boundary_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    boundary = gpd.read_file(boundary_path)

    if boundary.empty:
        raise ValueError(f"No boundary features found: {boundary_path}")

    if boundary.crs is None:
        raise ValueError(f"Boundary has no CRS: {boundary_path}")

    with rasterio.open(raster_path) as src:
        if src.crs is None:
            raise ValueError(f"Raster has no CRS: {raster_path}")

        raster = src.read(1, masked=True)
        extent = _raster_extent(src)
        raster_crs = src.crs

    boundary_plot = boundary.to_crs(raster_crs)

    valid_values = _get_valid_values(raster)

    if valid_values.size == 0:
        raise ValueError(f"Raster contains no valid pixels: {raster_path}")

    value_min = float(np.nanmin(valid_values))
    value_max = float(np.nanmax(valid_values))

    # Change maps show future minus baseline; zero is the neutral midpoint.
    is_change = _is_change_raster(raster_path)
    variable, variable_label, unit = _variable_info(raster_path)
    if is_change:
        colour_limit = max(abs(value_min), abs(value_max))
        plot_vmin, plot_vmax = -colour_limit, colour_limit
        plot_cmap = ("RdBu" if variable.lower() == "r20mm" else "RdBu_r") if cmap == DEFAULT_CMAP else cmap
    else:
        plot_vmin, plot_vmax = value_min, value_max
        plot_cmap = "RdYlBu" if variable.lower() == "r20mm" and cmap == DEFAULT_CMAP else cmap

    fig, ax = plt.subplots(figsize=(9, 8))

    image = ax.imshow(
        raster,
        extent=extent,
        origin="upper",
        cmap=plot_cmap,
        interpolation="nearest",
        vmin=plot_vmin,
        vmax=plot_vmax,
    )

    boundary_plot.plot(
        ax=ax,
        facecolor="none",
        edgecolor="black",
        linewidth=boundary_linewidth,
    )

    # Optional overlays: any number of user-provided shapefiles, each with its
    # own colour, transparency, outline width, and legend label. When no label
    # is supplied, use the shapefile name (without its extension).
    overlay_specs = _prepare_overlays(
        overlays,
        flood_extent_path=flood_extent_path,
        bushfire_prone_path=bushfire_prone_path,
        flood_color=flood_color,
        bushfire_color=bushfire_color,
        overlay_alpha=overlay_alpha,
    )
    overlay_legend = []
    for overlay in overlay_specs:
        overlay_legend.append(_overlay_layer(
            ax, overlay["path"], raster_crs,
            color=overlay["color"], label=overlay["label"],
            alpha=overlay["alpha"], linewidth=overlay["linewidth"],
        ))
    # Re-draw the study boundary on top of optional overlays.
    if overlay_legend:
        boundary_plot.plot(ax=ax, facecolor="none", edgecolor="black",
                           linewidth=boundary_linewidth, zorder=10)
        ax.legend(handles=overlay_legend, loc="upper right", fontsize=8)

    colour_bar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.035,
        pad=0.03,
    )
    if is_change:
        colour_bar.set_label(
            f"Change in {variable_label} ({_period_label(raster_path)} minus baseline; {unit})"
        )
    else:
        colour_bar.set_label(f"{variable_label} ({unit})")

    if title is None:
        title = _clean_raster_name(raster_path)

    ax.set_title(title, fontsize=12, pad=12)

    if raster_crs.to_epsg() == 4326:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    else:
        ax.set_xlabel("Easting")
        ax.set_ylabel("Northing")

    ax.set_aspect("equal")

    ax.text(
        0.01,
        0.01,
        f"Min: {value_min:.1f}   Max: {value_max:.1f}",
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="bottom",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "alpha": 0.75,
            "edgecolor": "none",
        },
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return output_path


# =============================================================================
# PLOT ALL ENSEMBLE RASTERS
# =============================================================================

def plot_ensemble_maps(
    *,
    input_root: str | Path,
    boundary_path: str | Path,
    climate_indices_folder: str = "Climate_Indicies",
    ensemble_folder_name: str = "Ensemble_Stats",
    output_folder_name: str = "Plot and Maps",
    variables: list[str] | None = None,
    scenarios: list[str] | None = None,
    statistics: list[str] | None = None,
    cmap: str = DEFAULT_CMAP,
    dpi: int = DEFAULT_DPI,
    overwrite: bool = False,
    show: bool = False,
    verbose: bool = True,
    flood_extent_path: str | Path | None = None,
    bushfire_prone_path: str | Path | None = None,
    flood_color: str = "lightskyblue",
    bushfire_color: str = "yellow",
    overlay_alpha: float = 0.22,
    overlays: list[dict] | None = None,
    create_spi12_changes: bool = True,
) -> list[Path]:
    """
    Create maps for ensemble-statistic GeoTIFFs produced by uncertainty.py.

    The function reads GeoTIFFs from:
        <input_root>/<climate_indices_folder>/<ensemble_folder_name>

    and saves PNG maps to:
        <input_root>/<climate_indices_folder>/<output_folder_name>

    Optional overlays accepts a list of dictionaries, e.g.:
        overlays=[
            {"path": "planning_areas.shp", "color": "yellow", "alpha": 0.15},
            {"path": "flood.shp", "label": "Flood extent", "color": "cyan"},
        ]
    If label is omitted, the shapefile's stem is used in the legend.
    The previous flood_extent_path/bushfire_prone_path arguments still work.
    """

    input_root = Path(input_root)
    boundary_path = Path(boundary_path)

    climate_indices_root = input_root / climate_indices_folder
    # Support the corrected folder name without breaking existing callers.
    if not climate_indices_root.exists() and climate_indices_folder == "Climate_Indicies":
        alternative_root = input_root / "Climate_Indices"
        if alternative_root.exists():
            climate_indices_root = alternative_root
    ensemble_folder = climate_indices_root / ensemble_folder_name
    output_folder = climate_indices_root / output_folder_name

    if not ensemble_folder.exists():
        raise FileNotFoundError(
            "Ensemble raster folder was not found:\n"
            f"{ensemble_folder}"
        )

    if not boundary_path.exists():
        raise FileNotFoundError(
            f"Boundary file was not found: {boundary_path}"
        )

    output_folder.mkdir(parents=True, exist_ok=True)

    # Generate missing SPI12 change rasters when both historical and future
    # ensemble statistics exist; leave existing uncertainty outputs unchanged.
    if create_spi12_changes:
        _create_missing_spi12_change_rasters(ensemble_folder, overwrite=False, verbose=verbose)

    # Discover both legacy flat outputs and new Min/Max/Mean/Change subfolders.
    raster_paths = sorted(ensemble_folder.rglob("*.tif"))

    if not raster_paths:
        raise FileNotFoundError(
            "No GeoTIFF files were found in:\n"
            f"{ensemble_folder}"
        )

    if verbose:
        print("=" * 100)
        print("NARCliM ENSEMBLE MAP WORKFLOW")
        print("=" * 100)
        print(f"[INPUT]  {ensemble_folder}")
        print(f"[OUTPUT] {output_folder}")
        print(f"[BOUNDARY] {boundary_path}")
        print(f"[RASTERS FOUND] {len(raster_paths)}")

    variable_filters = (
        {value.lower() for value in variables}
        if variables
        else None
    )

    scenario_filters = (
        {value.lower() for value in scenarios}
        if scenarios
        else None
    )

    statistic_filters = (
        {value.lower() for value in statistics}
        if statistics
        else None
    )

    saved_maps: list[Path] = []

    for raster_path in raster_paths:
        stem_lower = raster_path.stem.lower()

        if variable_filters:
            if not any(
                stem_lower.startswith(variable + "_")
                or stem_lower == variable
                for variable in variable_filters
            ):
                continue

        if scenario_filters:
            if not any(
                f"_{scenario}_" in stem_lower
                for scenario in scenario_filters
            ):
                continue

        if statistic_filters:

            statistic_matches = []

            for statistic in statistic_filters:

                statistic = statistic.lower()

                # New naming convention:
                #
                #   TXge35_Mean_ssp245_2040_2060.tif
                #
                new_pattern = (
                    f"_{statistic}_"
                )

                # Old naming convention:
                #
                #   TXge35_ssp245_mid_term_ID1_ensemble_mean.tif
                #
                old_pattern = (
                    f"_ensemble_{statistic}"
                )

                # Change files encode both the change type and its statistic.
                change_pattern = f"_change_{statistic}_"

                statistic_matches.append(
                    (change_pattern in stem_lower)
                    or (
                        new_pattern
                        in stem_lower
                    )
                    or stem_lower.endswith(
                        old_pattern
                    )
                )

            if not any(
                statistic_matches
            ):
                continue

        # Mirror the ensemble-statistic subfolder layout for PNG maps.
        category = raster_path.parent.name.capitalize()
        if category not in {"Min", "Max", "Mean", "Change"}:
            match = re.search(r"_Change_(Min|Mean|Max)_", raster_path.stem, re.IGNORECASE)
            if match:
                category = "Change"
            else:
                match = re.search(r"_(Min|Mean|Max)_", raster_path.stem, re.IGNORECASE)
                if not match:
                    match = re.search(r"_ensemble_(min|mean|max)$", raster_path.stem, re.IGNORECASE)
                category = match.group(1).capitalize() if match else "Other"

        output_path = output_folder / category / f"{raster_path.stem}.png"

        if output_path.exists() and not overwrite:
            if verbose:
                print(f"[SKIP] {output_path.name}")

            saved_maps.append(output_path)
            continue

        try:
            if verbose:
                print(f"[PLOT] {raster_path.name}")

            plot_ensemble_raster(
                raster_path=raster_path,
                boundary_path=boundary_path,
                output_path=output_path,
                cmap=cmap,
                dpi=dpi,
                show=show,
                flood_extent_path=flood_extent_path,
                bushfire_prone_path=bushfire_prone_path,
                flood_color=flood_color,
                bushfire_color=bushfire_color,
                overlay_alpha=overlay_alpha,
                overlays=overlays,
            )

            saved_maps.append(output_path)

            #if verbose:
                #print(f"[OK] {output_path.name}")

        except Exception as exc:
            print(
                "[ERROR] "
                f"{raster_path.name}: "
                f"{type(exc).__name__}: {exc}"
            )

    if verbose:
        print("\n" + "=" * 100)
        print("NARCliM ENSEMBLE MAP WORKFLOW FINISHED")
        print("=" * 100)
        print(f"[MAPS] {len(saved_maps)}")
        print(f"[OUTPUT] {output_folder}")

    return saved_maps