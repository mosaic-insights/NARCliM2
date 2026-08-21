from __future__ import annotations

from pathlib import Path
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio


# =============================================================================
# DEFAULT CONFIGURATION
# =============================================================================

DEFAULT_CMAP = "viridis"
DEFAULT_DPI = 300


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

    Old naming convention:
        TXge35_ssp245_mid_term_ID1_ensemble_mean.tif
        -> TXge35 | SSP245 | Mid Term | Ensemble Mean
    """

    stem = raster_path.stem

    # ==================================================================
    # NEW NAMING CONVENTION
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

    fig, ax = plt.subplots(figsize=(9, 8))

    image = ax.imshow(
        raster,
        extent=extent,
        origin="upper",
        cmap=cmap,
        interpolation="nearest",
        vmin=value_min,
        vmax=value_max,
    )

    boundary_plot.plot(
        ax=ax,
        facecolor="none",
        edgecolor="black",
        linewidth=boundary_linewidth,
    )

    colour_bar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.035,
        pad=0.03,
    )
    colour_bar.set_label("Climate index value")

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
) -> list[Path]:
    """
    Create maps for ensemble-statistic GeoTIFFs produced by uncertainty.py.

    The function reads GeoTIFFs from:
        <input_root>/<climate_indices_folder>/<ensemble_folder_name>

    and saves PNG maps to:
        <input_root>/<climate_indices_folder>/<output_folder_name>
    """

    input_root = Path(input_root)
    boundary_path = Path(boundary_path)

    climate_indices_root = input_root / climate_indices_folder
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

    raster_paths = sorted(ensemble_folder.glob("*.tif"))

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

                statistic_matches.append(
                    (
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

        output_path = output_folder / f"{raster_path.stem}.png"

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