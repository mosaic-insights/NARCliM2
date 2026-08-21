from __future__ import annotations

from pathlib import Path

import geopandas as gpd
from shapely.ops import unary_union

from .config import DOMAINS


DEFAULT_FOOTPRINT_PATH = (
    Path(__file__).resolve().parent
    / "Data"
    / "NARCliM2-0-SEAus-04_footprint.gpkg"
)


def load_seaus_footprint(
    footprint_path: str | Path | None = None,
):
    """Load the packaged NARCliM2-0-SEAus-04 footprint in EPSG:4326."""

    path = Path(footprint_path) if footprint_path else DEFAULT_FOOTPRINT_PATH

    if not path.exists():
        raise FileNotFoundError(
            "The SEAus domain footprint was not found:\n"
            f"{path}\n"
            "Run create_seaus_footprint.py once and save the GeoPackage "
            "inside src/narclim_workflow/Data/."
        )

    gdf = gpd.read_file(path)

    if gdf.empty:
        raise ValueError(f"The SEAus footprint contains no features: {path}")

    if gdf.crs is None:
        raise ValueError(f"The SEAus footprint has no CRS: {path}")

    gdf = gdf.to_crs("EPSG:4326")
    footprint = unary_union(gdf.geometry)

    if footprint.is_empty:
        raise ValueError(f"The SEAus footprint geometry is empty: {path}")

    return footprint


def select_domain_key(
    *,
    features,
    requested_domain_key: str,
    footprint_path: str | Path | None = None,
) -> tuple[str, str]:
    """
    Resolve auto, seaus_4km, or aus_18 to a configured domain key.

    Auto uses the 4 km domain only when every input feature is completely
    covered by the saved SEAus footprint. Otherwise it selects AUS-18.
    """

    valid_options = {"auto", *DOMAINS.keys()}

    if requested_domain_key not in valid_options:
        raise ValueError(
            f"Unknown domain_key '{requested_domain_key}'. "
            f"Available values: {sorted(valid_options)}"
        )

    if requested_domain_key != "auto":
        return requested_domain_key, (
            f"User explicitly selected '{requested_domain_key}'."
        )

    footprint = load_seaus_footprint(footprint_path)

    outside_ids = [
        str(feature.feature_id)
        for feature in features
        if not footprint.covers(feature.geometry)
    ]

    if not outside_ids:
        return "seaus_4km", (
            "All input features are completely covered by the "
            "NARCliM2-0-SEAus-04 footprint."
        )

    preview = ", ".join(outside_ids[:10])
    if len(outside_ids) > 10:
        preview += f", ... ({len(outside_ids)} features total)"

    return "aus_18", (
        "One or more input features extend outside the "
        "NARCliM2-0-SEAus-04 footprint. "
        f"Outside/intersecting feature IDs: {preview}"
    )
