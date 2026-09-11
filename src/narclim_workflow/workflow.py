from __future__ import annotations
import gc
from .domain import select_domain_key
from pathlib import Path
from typing import Iterable
import traceback

import pandas as pd
import numpy as np
import xarray as xr

from .catalog import build_catalog_url, read_catalog, select_files_for_window
from .config import DOMAINS, GCM_REALISATIONS, RCMS, SCENARIOS, VARIABLES
from .spatial import (
    choose_data_variable,
    extract_point,
    extract_polygon,
    read_features,
    subset_time,
)


def _log(message: str, *, verbose: bool) -> None:
    if verbose:
        print(message, flush=True)

def _validate_choices(*,domain_key,variables,gcms,scenarios,rcms,) -> None:
    # Allow automatic domain selection as well as explicit domains.
    valid_domain_options = {"auto", *DOMAINS.keys(),}
    if domain_key not in valid_domain_options:
        raise ValueError(f"Unknown domain_key {domain_key}. Available values: {sorted(valid_domain_options)}")
    unknown = set(variables) - set(VARIABLES)
    if unknown:
        raise ValueError(f"Unknown variables: {sorted(unknown)}")
    unknown = set(gcms) - set(GCM_REALISATIONS)
    if unknown:
        raise ValueError(f"Unknown GCMs: {sorted(unknown)}")
    unknown = set(scenarios) - set(SCENARIOS)
    if unknown:
        raise ValueError(f"Unknown scenarios: {sorted(unknown)}")
    unknown = set(rcms) - set(RCMS)
    if unknown:
        raise ValueError(f"Unknown RCMs: {sorted(unknown)}")

def _format_threshold_for_name(
    threshold_kmh: float,
) -> str:
    """
    Convert a numeric storm threshold to a safe variable-name component.

    Examples
    --------
    89.0 -> "89"
    70.5 -> "70p5"
    """

    value = float(threshold_kmh)

    if value <= 0:
        raise ValueError(
            "storm_threshold_kmh must be greater than zero."
        )

    if value.is_integer():
        return str(int(value))

    text = f"{value:.3f}".rstrip("0").rstrip(".")
    return text.replace(".", "p")


def _storm_index_name(
    threshold_kmh: float,
) -> str:
    """
    Return the derived storm-index variable/folder name.

    Example
    -------
    89.0 km/h -> StormDaysGT89
    70.5 km/h -> StormDaysGT70p5
    """

    return (
        "StormDaysGT"
        + _format_threshold_for_name(
            threshold_kmh
        )
    )

def _netcdf_safe_attrs(
    attrs: dict,
) -> dict:
    """
    Convert dataset attributes into values that can safely be written
    by the netCDF4 backend.

    NetCDF attributes do not support Python/NumPy Boolean values
    directly, so Booleans are converted to integers:

        True  -> 1
        False -> 0
    """

    safe_attrs = {}

    for key, value in attrs.items():

        if isinstance(
            value,
            (bool, np.bool_),
        ):
            safe_attrs[
                key
            ] = int(value)

        elif isinstance(
            value,
            Path,
        ):
            safe_attrs[
                key
            ] = str(value)

        else:
            safe_attrs[
                key
            ] = value

    return safe_attrs

def _derive_storm_days(
    ds: xr.Dataset,
    *,
    source_variable: str,
    threshold_kmh: float,
    output_variable_name: str,
) -> xr.Dataset:
    """
    Convert daily maximum near-surface wind speed into an annual count of
    days exceeding a user-defined wind-speed threshold.

    Source
    ------
    NARCliM sfcWindmax:
        Daily Maximum Near-Surface Wind Speed
        units = m s-1

    Derived index
    -------------
    Number of days in each year where daily maximum near-surface wind speed
    is greater than ``threshold_kmh``.

    Notes
    -----
    The threshold is supplied in km/h and converted internally to m/s.

    The calculation uses the actual time values in the source NetCDF and
    therefore does not assume that every climate-model year contains
    365 days.
    """

    if ds is None:
        raise ValueError(
            "_derive_storm_days received ds=None. "
            "The spatial extraction step did not return a dataset."
        )

    if not isinstance(
        ds,
        xr.Dataset,
    ):
        raise TypeError(
            "_derive_storm_days expected an xarray.Dataset, "
            f"but received {type(ds)}."
        )

    if source_variable is None:
        raise ValueError(
            "source_variable is None."
        )

    if source_variable not in ds.data_vars:
        raise KeyError(
            f"Source wind variable '{source_variable}' "
            "was not found in the dataset. "
            f"Available data variables: {list(ds.data_vars)}"
        )

    if threshold_kmh is None:
        raise ValueError(
            "threshold_kmh is None."
        )

    threshold_kmh = float(
        threshold_kmh
    )

    if threshold_kmh <= 0:
        raise ValueError(
            "threshold_kmh must be greater than zero."
        )

    if "time" not in ds[source_variable].dims:
        raise ValueError(
            f"'{source_variable}' does not contain a time dimension."
        )

    wind = ds[
        source_variable
    ]

    # --------------------------------------------------------------
    # Validate units.
    # --------------------------------------------------------------

    units = str(
        wind.attrs.get(
            "units",
            "",
        )
    ).strip()

    accepted_units = {
        "m s-1",
        "m s^-1",
        "m/s",
        "m s**-1",
    }

    if units not in accepted_units:
        raise ValueError(
            "Storm-day derivation expects sfcWindmax in metres per second. "
            f"Found units={units!r}."
        )

    # --------------------------------------------------------------
    # Convert threshold from km/h to m/s.
    # --------------------------------------------------------------

    threshold_ms = (
        threshold_kmh
        / 3.6
    )

    # --------------------------------------------------------------
    # Preserve the original valid-data footprint.
    #
    # NaN > threshold evaluates to False, so without this mask an
    # entirely missing pixel could incorrectly become zero storm days.
    # --------------------------------------------------------------

    valid_days = (
        wind.notnull()
        .sum(
            dim="time"
        )
    )

    exceedance_days = (
        wind
        > threshold_ms
    )

    annual_count = (
        exceedance_days
        .sum(
            dim="time"
        )
        .where(
            valid_days > 0
        )
        .astype(
            "float32"
        )
    )

    # --------------------------------------------------------------
    # Each NCI daily sfcWindmax source file represents one model year.
    # Keep a one-element time dimension so uncertainty.py can later
    # concatenate years and calculate the temporal average.
    # --------------------------------------------------------------

    first_time = (
        ds["time"]
        .values[0]
    )

    annual_count = (
        annual_count
        .expand_dims(
            time=[
                first_time
            ]
        )
    )

    annual_count.name = (
        output_variable_name
    )

    annual_count.attrs = {
        "long_name": (
            "Annual number of days with daily maximum near-surface "
            f"wind speed greater than {threshold_kmh:g} km/h"
        ),
        "description": (
            "Count of days where NARCliM daily maximum near-surface "
            f"wind speed (sfcWindmax) exceeds {threshold_kmh:g} km/h "
            f"({threshold_ms:.6f} m/s)."
        ),
        # A count is dimensionless. Keeping units='1' avoids xarray
        # interpreting the variable as a timedelta.
        "units": "1",
        "count_units": "days",
        "temporal_basis": "annual",
        "source_variable": source_variable,
        "source_variable_units": units,
        "threshold_kmh": threshold_kmh,
        "threshold_ms": float(threshold_ms),
        "cell_methods": "time: maximum within days time: sum over days",
        "derived_index": output_variable_name,
    }

    # --------------------------------------------------------------
    # Build output dataset.
    # --------------------------------------------------------------

    derived_ds = (
        annual_count
        .to_dataset()
    )

    # Preserve useful spatial/scalar coordinates from the extracted
    # NARCliM dataset.
    for coordinate_name in (
        "rlon",
        "rlat",
        "lon",
        "lat",
        "height",
        "crs",
    ):

        if (
            coordinate_name
            in ds.variables
            and coordinate_name
            not in derived_ds.variables
        ):

            derived_ds[
                coordinate_name
            ] = ds[
                coordinate_name
            ]

            if coordinate_name in ds.coords:
                derived_ds = (
                    derived_ds
                    .set_coords(
                        coordinate_name
                    )
                )

    # Preserve global NARCliM metadata.
    derived_ds.attrs = (
        ds.attrs.copy()
    )

    derived_ds.attrs.update(
        {
            "variable_id": output_variable_name,
            "derived_from": source_variable,
            "storm_threshold_kmh": threshold_kmh,
            "storm_threshold_ms": float(threshold_ms),
            "derived_frequency": "yr",
        }
    )

    return derived_ds


def run_workflow(
    *,
    vector_path,
    output_root,
    time_windows,
    variables,
    gcms,
    scenarios,
    rcms,
    domain_key="seaus_4km",
    id_field=None,
    overwrite=False,
    engine="netcdf4",
    verbose=True,
    print_traceback=False,
    storm_threshold_kmh: float | None = None,
    save_storm_source: bool = True,
) -> pd.DataFrame:
    """Run the NARCliM workflow with detailed progress messages."""

    _validate_choices(
        domain_key=domain_key,
        variables=variables,
        gcms=gcms,
        scenarios=scenarios,
        rcms=rcms,
    )

    output_root = Path(output_root) / "NARCliM Data"
    output_root.mkdir(parents=True, exist_ok=True)

    _log("=" * 100, verbose=verbose)
    _log("NARCliM WORKFLOW STARTED", verbose=verbose)
    _log("=" * 100, verbose=verbose)
    _log(f"[CONFIG] Vector file : {vector_path}", verbose=verbose)
    _log(f"[CONFIG] Output root : {output_root}", verbose=verbose)
    _log(f"[CONFIG] Domain mode : {domain_key}",verbose=verbose,)
    _log(f"[CONFIG] Variables   : {variables}", verbose=verbose)
    _log(f"[CONFIG] Scenarios   : {scenarios}", verbose=verbose)
    _log(f"[CONFIG] GCMs        : {gcms}", verbose=verbose)
    _log(f"[CONFIG] RCMs        : {rcms}", verbose=verbose)
    _log(f"[CONFIG] Time windows: {time_windows}", verbose=verbose)
    _log(f"[CONFIG] xarray engine: {engine}", verbose=verbose)
    _log(
        f"[CONFIG] Storm threshold: "
        f"{storm_threshold_kmh if storm_threshold_kmh is not None else 'config default'} km/h",
        verbose=verbose,
    )
    _log(
        f"[CONFIG] Save raw storm source (sfcWindmax): {save_storm_source}",
        verbose=verbose,
    )

    _log("\n[STEP 1] Reading vector features...", verbose=verbose)
    features = read_features(vector_path, id_field=id_field)
    selected_domain_key, domain_reason = select_domain_key(features=features,requested_domain_key=domain_key,)
    domain = DOMAINS[selected_domain_key]
    _log(f"[DOMAIN REQUEST] {domain_key}",verbose=verbose,)

    _log(
        f"[DOMAIN SELECTED] "
        f"{selected_domain_key} -> {domain}",
        verbose=verbose,)

    _log(
        f"[DOMAIN REASON] {domain_reason}",
        verbose=verbose,)    
    _log(f"[OK] Loaded {len(features)} feature(s).", verbose=verbose)
    for feature in features:
        _log(
            f"     - ID={feature.feature_id}; geometry={feature.geometry.geom_type}; "
            f"bounds={tuple(round(v, 6) for v in feature.geometry.bounds)}",
            verbose=verbose,
        )

    records = []
    combination_number = 0

    for variable_key in variables:
        spec = VARIABLES[variable_key]

        _log(
            "\n" + "#" * 100,
            verbose=verbose,
        )

        _log(
            f"[VARIABLE] {variable_key}: {spec.description}",
            verbose=verbose,
        )
        for scenario in scenarios:
            if scenario not in spec.available_scenarios:
                message = (
                    f"{variable_key} is not available for scenario "
                    f"'{scenario}'. Continuing to the next scenario.")
                _log(
                    f"[DATA NOT AVAILABLE] {message}",
                    verbose=verbose,)
                records.append(
                    {
                        "status": "data_not_available",
                        "variable": variable_key,
                        "scenario": scenario,
                        "window": None,
                        "gcm": None,
                        "rcm": None,
                        "catalog_url": None,
                        "source_url": None,
                        "feature_id": None,
                        "output": None,
                        "message": message,
                    }
                )

                continue
            for window_name, (start_year, end_year) in time_windows.items():
                # Historical NARCliM data end in 2014.
                if scenario == "historical" and start_year > 2014:
                    message = (
                        f"Time window '{window_name}' ({start_year}-{end_year}) "
                        f"does not overlap the historical period.")
                    _log(
                        f"[TIME WINDOW NOT AVAILABLE] {message}",
                        verbose=verbose,)
                    records.append(
                        {
                            "status": "time_window_not_available",
                            "variable": variable_key,
                            "scenario": scenario,
                            "window": window_name,
                            "gcm": None,
                            "rcm": None,
                            "catalog_url": None,
                            "source_url": None,
                            "feature_id": None,
                            "output": None,
                            "message": message,
                        }
                    )
                    continue
                # NARCliM SSP projections begin in 2015.
                if scenario != "historical" and end_year < 2015:
                    message = (
                        f"Time window '{window_name}' ({start_year}-{end_year}) "
                        f"does not overlap projection scenario '{scenario}'."
                    )
                    _log(
                        f"[TIME WINDOW NOT AVAILABLE] {message}",
                        verbose=verbose,
                    )
                    records.append(
                        {
                            "status": "time_window_not_available",
                            "variable": variable_key,
                            "scenario": scenario,
                            "window": window_name,
                            "gcm": None,
                            "rcm": None,
                            "catalog_url": None,
                            "source_url": None,
                            "feature_id": None,
                            "output": None,
                            "message": message,
                        }
                    )
                    continue

                for gcm in gcms:
                    for rcm in rcms:
                        combination_number += 1
                        _log("\n" + "-" * 100, verbose=verbose)
                        _log(
                            f"[COMBINATION {combination_number}] variable={variable_key} | "
                            f"scenario={scenario} | window={window_name} ({start_year}-{end_year}) | "
                            f"GCM={gcm} | RCM={rcm}",
                            verbose=verbose,
                        )

                        catalog_url = build_catalog_url(
                            domain=domain,
                            gcm=gcm,
                            scenario=scenario,
                            rcm=rcm,
                            variable=spec,
                        )
                        #_log(f"[CATALOG URL] {catalog_url}", verbose=verbose)

                        try:
                            #_log("[CATALOG] Reading catalog.xml...", verbose=verbose)
                            all_files = read_catalog(catalog_url)
                            _log(f"[CATALOG] Found {len(all_files)} NetCDF file(s).", verbose=verbose)
                            """
                            for item in all_files[:3]:
                                _log(f"          {item.name}", verbose=verbose)
                            if len(all_files) > 3:
                                _log(f"          ... and {len(all_files)-3} more", verbose=verbose)
                            """
                            files = select_files_for_window(all_files,start_year=start_year,end_year=end_year,)
                            _log(
                                f"[TIME FILTER] {len(files)} file(s) overlap {start_year}-{end_year}.",
                                verbose=verbose,
                            )

                            if not files:
                                message = "No source files overlap the requested time window."
                                _log(f"[WARNING] {message}", verbose=verbose)
                                records.append({
                                    "status": "no_matching_files",
                                    "variable": variable_key,
                                    "scenario": scenario,
                                    "window": window_name,
                                    "gcm": gcm,
                                    "rcm": rcm,
                                    "catalog_url": catalog_url,
                                    "source_url": None,
                                    "feature_id": None,
                                    "output": None,
                                    "message": message,
                                })
                                continue
                        except Exception as exc:
                            error_text = str(exc)
                            error_type = type(exc).__name__
                            # Missing catalogue/data should not stop the complete workflow.
                            likely_unavailable = (
                                "404" in error_text
                                or "No NetCDF files found" in error_text
                                or "not found" in error_text.lower()
                            )
                            if likely_unavailable:
                                status = "data_not_available"
                                message = (
                                    f"No data were found for {variable_key}, {scenario}, "
                                    f"{window_name}, {gcm}, {rcm}. "
                                    f"Continuing to the next combination. "
                                    f"Details: {error_type}: {error_text}"
                                )
                                _log(
                                    f"[DATA NOT AVAILABLE] {message}",
                                    verbose=verbose,
                                )
                            else:
                                status = "catalog_error"
                                message = f"{error_type}: {error_text}"
                                _log(
                                    f"[CATALOG ERROR] {message}",
                                    verbose=verbose,
                                )
                                if print_traceback:
                                    traceback.print_exc()
                            records.append(
                                {
                                    "status": status,
                                    "variable": variable_key,
                                    "scenario": scenario,
                                    "window": window_name,
                                    "gcm": gcm,
                                    "rcm": rcm,
                                    "catalog_url": catalog_url,
                                    "source_url": None,
                                    "feature_id": None,
                                    "output": None,
                                    "message": message,
                                }
                            )
                            continue
                        for file_index, source_file in enumerate(files, start=1):
                            _log(
                                f"\n[SOURCE {file_index}/{len(files)}] {source_file.name}",
                                verbose=verbose,
                            )
                            #_log(f"[OPENDAP URL] {source_file.opendap_url}", verbose=verbose)

                            try:
                                #_log("[OPEN] Opening remote dataset...", verbose=verbose)
                                with xr.open_dataset(
                                    source_file.opendap_url,
                                    engine=engine,
                                    decode_times=True,

                                    # Do not create a Dask graph for these small spatial subsets.
                                    chunks=None,

                                    # Avoid retaining remote arrays in xarray's in-memory cache.
                                    cache=False,
                                ) as source_ds:
                                    #_log("[OPEN] Dataset opened.", verbose=verbose)
                                    #_log(f"[DATASET] Sizes: {dict(source_ds.sizes)}", verbose=verbose)
                                    #_log(f"[DATASET] Data variables: {list(source_ds.data_vars)}", verbose=verbose)
                                    #_log(f"[DATASET] Coordinates: {list(source_ds.coords)}", verbose=verbose)

                                    data_var = choose_data_variable(source_ds, spec.name)
                                    #_log(f"[VARIABLE] Using '{data_var}'.", verbose=verbose)

                                    temporal_ds = subset_time(source_ds, start_year, end_year)
                                    #_log(f"[TIME SUBSET] Sizes: {dict(temporal_ds.sizes)}",verbose=verbose,)

                                    if any(size == 0 for size in temporal_ds.sizes.values()):
                                        message = "Subset contains a zero-length dimension."
                                        _log(f"[WARNING] {message}", verbose=verbose)
                                        records.append({
                                            "status": "empty_subset",
                                            "variable": variable_key,
                                            "scenario": scenario,
                                            "window": window_name,
                                            "gcm": gcm,
                                            "rcm": rcm,
                                            "catalog_url": catalog_url,
                                            "source_url": source_file.opendap_url,
                                            "feature_id": None,
                                            "output": None,
                                            "message": message,
                                        })
                                        continue

                                    for feature_index, feature in enumerate(
                                        features,
                                        start=1,
                                    ):

                                        geom = feature.geometry

                                        out_ds = None

                                        _log(
                                            f"[FEATURE {feature_index}/{len(features)}] "
                                            f"ID={feature.feature_id}; "
                                            f"type={geom.geom_type}",
                                            verbose=verbose,
                                        )

                                        try:

                                            # ==============================================================
                                            # SPATIAL EXTRACTION
                                            # ==============================================================

                                            if geom.geom_type == "Point":

                                                #_log("[SPATIAL] Extracting nearest grid cell...",verbose=verbose,)

                                                out_ds = extract_point(
                                                    temporal_ds,
                                                    geom,
                                                    data_var,
                                                )

                                            elif geom.geom_type in {
                                                "Polygon",
                                                "MultiPolygon",
                                            }:

                                                #_log("[SPATIAL] Applying polygon subset/mask...",verbose=verbose,)

                                                out_ds = extract_polygon(
                                                    temporal_ds,
                                                    geom,
                                                    data_var,
                                                    verbose=False,
                                                )

                                            else:

                                                raise ValueError(
                                                    f"Unsupported geometry type: {geom.geom_type}"
                                                )


                                            # ==============================================================
                                            # VALIDATE SPATIAL OUTPUT
                                            # ==============================================================

                                            if out_ds is None:

                                                raise RuntimeError(
                                                    "Spatial extraction returned None. "
                                                    "Expected an xarray.Dataset."
                                                )

                                            if not isinstance(
                                                out_ds,
                                                xr.Dataset,
                                            ):

                                                raise TypeError(
                                                    "Spatial extraction returned an unexpected type: "
                                                    f"{type(out_ds)}"
                                                )


                                            #_log(f"[SPATIAL] Output sizes: {dict(out_ds.sizes)}",verbose=verbose,)


                                            # ==============================================================
                                            # DERIVED CLIMATE INDICES
                                            # ==============================================================

                                            # By default, outputs retain the user-requested variable key.
                                            output_variable_key = variable_key

                                            is_storm_request = (
                                                spec.name == "sfcWindmax"
                                                and getattr(
                                                    spec,
                                                    "derived_variable",
                                                    None,
                                                )
                                                is not None
                                            )

                                            if is_storm_request:

                                                # ----------------------------------------------------------
                                                # Resolve the user-defined threshold.
                                                #
                                                # Priority:
                                                #   1. run_workflow(storm_threshold_kmh=...)
                                                #   2. VariableSpec.threshold_kmh default
                                                # ----------------------------------------------------------

                                                resolved_storm_threshold = (
                                                    storm_threshold_kmh
                                                    if storm_threshold_kmh is not None
                                                    else getattr(
                                                        spec,
                                                        "threshold_kmh",
                                                        None,
                                                    )
                                                )

                                                if resolved_storm_threshold is None:
                                                    raise ValueError(
                                                        "Storm processing requires a threshold. "
                                                        "Set storm_threshold_kmh in run_workflow() "
                                                        "or threshold_kmh in VariableSpec."
                                                    )

                                                resolved_storm_threshold = float(
                                                    resolved_storm_threshold
                                                )

                                                derived_storm_name = (
                                                    _storm_index_name(
                                                        resolved_storm_threshold
                                                    )
                                                )

                                                # ----------------------------------------------------------
                                                # Load the small ACT/study-area subset once.
                                                #
                                                # This source subset is then:
                                                #   1. optionally saved as raw sfcWindmax;
                                                #   2. used to derive annual storm-day counts.
                                                # ----------------------------------------------------------

                                                out_ds.load()

                                                rcm_short = {
                                                    "NARCliM2-0-WRF412R3": "R3",
                                                    "NARCliM2-0-WRF412R5": "R5",
                                                }

                                                year_part = (
                                                    Path(
                                                        source_file.name
                                                    )
                                                    .stem
                                                    .split("_")[-1]
                                                )

                                                # ----------------------------------------------------------
                                                # SAVE ORIGINAL DAILY STORM SOURCE
                                                #
                                                # This makes sfcWindmax behave like the other downloaded
                                                # source variables (e.g. precipitation and temperature):
                                                #
                                                # NARCliM Data/
                                                #   sfcWindmax/
                                                #     historical/
                                                #       baseline/
                                                #         GCM/R3/ID1/
                                                # ----------------------------------------------------------

                                                if save_storm_source:

                                                    source_folder = (
                                                        output_root
                                                        / "sfcWindmax"
                                                        / scenario
                                                        / window_name
                                                        / gcm
                                                        / rcm_short.get(
                                                            rcm,
                                                            rcm,
                                                        )
                                                        / str(
                                                            feature.feature_id
                                                        )
                                                    )

                                                    source_folder.mkdir(
                                                        parents=True,
                                                        exist_ok=True,
                                                    )

                                                    source_output_path = (
                                                        source_folder
                                                        / (
                                                            f"sfcWindmax_"
                                                            f"{year_part}.nc"
                                                        )
                                                    )

                                                    if (
                                                        source_output_path.exists()
                                                        and not overwrite
                                                    ):
                                                        source_status = (
                                                            "skipped_existing"
                                                        )
                                                    else:

                                                        source_to_save = (
                                                            out_ds.copy(
                                                                deep=False
                                                            )
                                                        )

                                                        source_to_save.attrs.update(
                                                            _netcdf_safe_attrs(
                                                                {
                                                                    "source_catalog":
                                                                        catalog_url,

                                                                    "source_opendap":
                                                                        source_file.opendap_url,

                                                                    "requested_start_year":
                                                                        start_year,

                                                                    "requested_end_year":
                                                                        end_year,

                                                                    "feature_id":
                                                                        feature.feature_id,

                                                                    "saved_as_original_storm_source":
                                                                        True,
                                                                }
                                                            )
                                                        )

                                                        source_temp_path = (
                                                            source_output_path
                                                            .with_suffix(
                                                                ".tmp.nc"
                                                            )
                                                        )

                                                        source_temp_path.unlink(
                                                            missing_ok=True
                                                        )

                                                        source_to_save.to_netcdf(
                                                            source_temp_path,
                                                            engine="netcdf4",
                                                            mode="w",
                                                        )

                                                        source_temp_path.replace(
                                                            source_output_path
                                                        )

                                                        source_status = "saved"

                                                    records.append(
                                                        {
                                                            "status":
                                                                source_status,
                                                            "variable":
                                                                "sfcWindmax",
                                                            "scenario":
                                                                scenario,
                                                            "window":
                                                                window_name,
                                                            "gcm":
                                                                gcm,
                                                            "rcm":
                                                                rcm,
                                                            "catalog_url":
                                                                catalog_url,
                                                            "source_url":
                                                                source_file.opendap_url,
                                                            "feature_id":
                                                                feature.feature_id,
                                                            "output":
                                                                str(
                                                                    source_output_path
                                                                ),
                                                            "message":
                                                                (
                                                                    "Original daily "
                                                                    "storm source subset."
                                                                ),
                                                        }
                                                    )

                                                    #_log("[STORM SOURCE SAVED] " f"{source_output_path}", verbose=verbose,)

                                                # ----------------------------------------------------------
                                                # DERIVE ANNUAL STORM-DAY INDEX
                                                # ----------------------------------------------------------

                                                #_log("[DERIVED INDEX] Calculating annual number of days with daily maximum wind speed > "
                                                    #f"{resolved_storm_threshold:g} km/h " f"as {derived_storm_name}...",verbose=verbose,)

                                                out_ds = _derive_storm_days(
                                                    out_ds,
                                                    source_variable=data_var,
                                                    threshold_kmh=(
                                                        resolved_storm_threshold
                                                    ),
                                                    output_variable_name=(
                                                        derived_storm_name
                                                    ),
                                                )

                                                data_var = (
                                                    derived_storm_name
                                                )

                                                output_variable_key = (
                                                    derived_storm_name
                                                )


                                            # ==============================================================
                                            # SAFETY CHECK FOR POLYGON EXTRACTION
                                            # ==============================================================

                                            if geom.geom_type in {
                                                "Polygon",
                                                "MultiPolygon",
                                            }:

                                                source_rows = (
                                                    temporal_ds.sizes.get(
                                                        "rlat"
                                                    )
                                                )

                                                source_columns = (
                                                    temporal_ds.sizes.get(
                                                        "rlon"
                                                    )
                                                )

                                                output_rows = (
                                                    out_ds.sizes.get(
                                                        "rlat"
                                                    )
                                                )

                                                output_columns = (
                                                    out_ds.sizes.get(
                                                        "rlon"
                                                    )
                                                )

                                                full_grid_retained = (
                                                    source_rows is not None
                                                    and source_columns is not None
                                                    and output_rows == source_rows
                                                    and output_columns == source_columns
                                                )

                                                if full_grid_retained:

                                                    raise RuntimeError(
                                                        "Polygon extraction retained the complete "
                                                        f"{source_rows} × {source_columns} grid. "
                                                        "This probably indicates a failed "
                                                        "curvilinear-grid subset."
                                                    )

                                            rcm_short = {
                                                "NARCliM2-0-WRF412R3": "R3",
                                                "NARCliM2-0-WRF412R5": "R5",
                                            }
                                            folder = (
                                                output_root
                                                / output_variable_key
                                                / scenario
                                                / window_name
                                                / gcm
                                                / rcm_short.get(rcm, rcm)
                                                / str(feature.feature_id)
                                            )
                                            folder.mkdir(
                                                parents=True,
                                                exist_ok=True,
                                            )
                                            year_part = (
                                                Path(source_file.name)
                                                .stem
                                                .split("_")[-1]
                                            )
                                            filename = (
                                                f"{output_variable_key}_"
                                                f"{year_part}.nc"
                                            )
                                            output_path = folder / filename
                                            _log(f"[OUTPUT SAVED] {output_path}",verbose=verbose,)
                                            if output_path.exists() and not overwrite:
                                                status = "skipped_existing"
                                                #_log("[SKIP] Output already exists.",verbose=verbose,)
                                            else:
                                                out_ds.attrs.update(
                                                    {
                                                        "source_catalog": catalog_url,
                                                        "source_opendap": (source_file.opendap_url),
                                                        "requested_start_year": (start_year),
                                                        "requested_end_year": (end_year),
                                                        "feature_id": (feature.feature_id),
                                                    })
                                                #_log("[LOAD] Loading remote subset ""into memory...",verbose=verbose,)
                                                # The safety check above ensures this is a small subset.
                                                out_ds.load()
                                                #_log("[LOAD] Remote subset loaded.",verbose=verbose,)
                                                #_log("[SAVE] Writing NetCDF...",verbose=verbose,)
                                                temporary_path = (output_path.with_suffix(".tmp.nc"))
                                                temporary_path.unlink(missing_ok=True)
                                                out_ds.to_netcdf(temporary_path,engine="netcdf4",mode="w",)
                                                temporary_path.replace(output_path)
                                                status = "saved"
                                                size_mb = (output_path.stat().st_size / (1024 ** 2))
                                                #_log(f"[OK] Saved {size_mb:.2f} MB.",verbose=verbose,)
                                            records.append(
                                                {
                                                    "status": status,
                                                    "variable": output_variable_key,
                                                    "scenario": scenario,
                                                    "window": window_name,
                                                    "gcm": gcm,
                                                    "rcm": rcm,
                                                    "catalog_url": catalog_url,
                                                    "source_url": (source_file.opendap_url),
                                                    "feature_id": (feature.feature_id),
                                                    "output": str(output_path),
                                                    "message": "",
                                                })
                                        except Exception as exc:
                                            message = (f"{type(exc).__name__}: {exc}")
                                            _log(f"[PROCESSING ERROR] {message}",verbose=verbose,)
                                            if print_traceback:
                                                traceback.print_exc()
                                            records.append(
                                                {
                                                    "status": "processing_error",
                                                    "variable": variable_key,
                                                    "scenario": scenario,
                                                    "window": window_name,
                                                    "gcm": gcm,
                                                    "rcm": rcm,
                                                    "catalog_url": catalog_url,
                                                    "source_url": (source_file.opendap_url),
                                                    "feature_id": (feature.feature_id),
                                                    "output": None,
                                                    "message": message,
                                                })
                                        finally:
                                            # Explicitly release the derived Dataset. The source Dataset is
                                            # closed automatically by the enclosing `with` statement.
                                            if out_ds is not None:
                                                try:
                                                    out_ds.close()
                                                except Exception:
                                                    pass

                                            del out_ds
                                            gc.collect()

                            except Exception as exc:
                                # Catch errors that occur while opening or processing
                                # the current remote source file.
                                message = f"{type(exc).__name__}: {exc}"

                                _log(
                                    f"[SOURCE PROCESSING ERROR] {message}",
                                    verbose=verbose,
                                )

                                if print_traceback:
                                    traceback.print_exc()

                                records.append(
                                    {
                                        "status": "processing_error",
                                        "variable": variable_key,
                                        "scenario": scenario,
                                        "window": window_name,
                                        "gcm": gcm,
                                        "rcm": rcm,
                                        "catalog_url": catalog_url,
                                        "source_url": source_file.opendap_url,
                                        "feature_id": None,
                                        "output": None,
                                        "message": message,
                                    }
                                )

                            finally:
                                # source_ds is closed by the with block.
                                gc.collect()

    manifest = pd.DataFrame(records)
    manifest_path = output_root / "narclim_download_manifest.csv"
    manifest.to_csv(manifest_path,index=False,)
    _log(f"[MANIFEST] Saved: {manifest_path}",verbose=verbose,)

    _log("\n" + "=" * 100, verbose=verbose)
    _log("NARCliM WORKFLOW FINISHED", verbose=verbose)
    _log("=" * 100, verbose=verbose)
    _log(f"[MANIFEST] {manifest_path}", verbose=verbose)

    if not manifest.empty:
        _log("\n[STATUS SUMMARY]", verbose=verbose)
        _log(manifest["status"].value_counts(dropna=False).to_string(), verbose=verbose)

        errors = manifest.loc[manifest["status"].isin(
            ["catalog_error", "processing_error", "no_matching_files", "empty_subset"]
        )]
        if not errors.empty:
            _log("\n[ERROR/WARNING PREVIEW]", verbose=verbose)
            for _, row in errors.head(10).iterrows():
                _log(
                    f"  - {row['status']} | {row.get('variable')} | "
                    f"{row.get('scenario')} | {row.get('gcm')} | "
                    f"{row.get('rcm')} | {row.get('message')}",
                    verbose=verbose,
                )

    return manifest