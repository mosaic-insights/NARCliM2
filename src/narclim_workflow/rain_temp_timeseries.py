from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

SOURCE_VARIABLES = ("prAdjust", "tasmaxAdjust", "tasminAdjust")
ALL_VARIABLES = ("prAdjust", "tasmaxAdjust", "tasminAdjust", "tasmeanAdjust")
FREQUENCIES = ("Monthly", "Seasonal", "Annual")
STATISTICS = ("min", "mean", "max")

SEASON_ORDER = {"DJF": 1, "MAM": 2, "JJA": 3, "SON": 4}
SEASON_REPRESENTATIVE_MONTH = {"DJF": 2, "MAM": 5, "JJA": 8, "SON": 11}

VARIABLE_LABELS = {
    "prAdjust": "Rainfall",
    "tasmaxAdjust": "Maximum temperature",
    "tasminAdjust": "Minimum temperature",
    "tasmeanAdjust": "Mean temperature",
}

VARIABLE_UNITS = {
    "prAdjust": "mm",
    "tasmaxAdjust": "degC",
    "tasminAdjust": "degC",
    "tasmeanAdjust": "degC",
}


def _log(message: str, *, verbose: bool) -> None:
    if verbose:
        print(message, flush=True)


def _normalise_root(input_root: str | Path) -> tuple[Path, Path]:
    path = Path(input_root)
    if path.name.lower() == "narclim data":
        data_root = path
        project_root = path.parent
    else:
        project_root = path
        data_root = path / "NARCliM Data"

    if not data_root.exists():
        direct = [project_root / v for v in SOURCE_VARIABLES]
        if any(p.exists() for p in direct):
            data_root = project_root
        else:
            raise FileNotFoundError(
                "Could not find NARCliM Data or variable folders under input_root."
            )
    return project_root, data_root


def _safe_name(value: object) -> str:
    text = str(value).strip()
    for old, new in [(" ", "_"), ("/", "-"), ("\\", "-"), (":", "-")]:
        text = text.replace(old, new)
    return text


def _choose_variable(ds: xr.Dataset, expected: str) -> str:
    if expected in ds.data_vars:
        return expected
    excluded = {"crs", "time_bnds", "lat", "lon", "rlat", "rlon", "height"}
    candidates = [
        name for name, da in ds.data_vars.items()
        if name not in excluded and "time" in da.dims
    ]
    if len(candidates) == 1:
        return candidates[0]
    raise KeyError(
        f"Could not identify '{expected}'. Available variables: {list(ds.data_vars)}"
    )


def _spatial_dims(da: xr.DataArray) -> list[str]:
    excluded = {"time", "ensemble_member", "season_year", "season", "_group"}
    return [dim for dim in da.dims if dim not in excluded]


def _temperature_to_celsius(da: xr.DataArray, *, variable: str) -> xr.DataArray:
    units = str(da.attrs.get("units", "")).strip().lower()
    kelvin_units = {"k", "kelvin", "degrees_k", "degree_k"}
    celsius_units = {
        "c", "°c", "degc", "degree_celsius", "degrees_celsius",
        "degree celsius", "degrees celsius",
    }
    if units in kelvin_units:
        out = da - 273.15
        conversion = "converted from Kelvin to degC"
    elif units in celsius_units:
        out = da.copy(deep=False)
        conversion = "source already in degC"
    else:
        raise ValueError(
            f"{variable}: unsupported temperature units {units!r}."
        )
    out.attrs = da.attrs.copy()
    out.attrs["units"] = "degC"
    out.attrs["unit_conversion"] = conversion
    return out


def _precipitation_to_mm_per_day(da: xr.DataArray) -> xr.DataArray:
    units_raw = str(da.attrs.get("units", "")).strip()
    units = units_raw.lower().replace("**", "^")
    flux_units = {
        "kg m-2 s-1", "kg m^-2 s^-1", "kg/m2/s",
        "mm s-1", "mm s^-1", "mm/s",
    }
    daily_units = {
        "mm", "mm d-1", "mm day-1", "mm/day",
        "mm d^-1", "mm day^-1", "kg m-2 d-1", "kg m^-2 d^-1",
    }
    if units in flux_units:
        out = da * 86400.0
        conversion = "daily flux multiplied by 86400 s/day"
    elif units in daily_units:
        out = da.copy(deep=False)
        conversion = "source interpreted as daily precipitation amount"
    else:
        raise ValueError(
            f"prAdjust: unsupported precipitation units {units_raw!r}."
        )
    out.attrs = da.attrs.copy()
    out.attrs["units"] = "mm/day"
    out.attrs["unit_conversion"] = conversion
    return out


@dataclass(frozen=True)
class MemberKey:
    gcm: str
    rcm: str
    feature_id: str

    @property
    def label(self) -> str:
        return f"{self.gcm}__{self.rcm}"


def _discover_member_files(
    *, data_root: Path, variable: str, scenario: str, window: str,
    gcms: Iterable[str] | None = None,
    rcms: Iterable[str] | None = None,
    feature_ids: Iterable[str] | None = None,
) -> dict[MemberKey, list[Path]]:
    base = data_root / variable / scenario / window
    if not base.exists():
        return {}

    allowed_gcms = None if gcms is None else set(gcms)
    allowed_features = None if feature_ids is None else {str(v) for v in feature_ids}
    requested_aliases = None
    if rcms is not None:
        requested_aliases = set()
        for item in rcms:
            text = str(item)
            requested_aliases.add(text)
            requested_aliases.add(text.replace("NARCliM2-0-WRF412", ""))

    result: dict[MemberKey, list[Path]] = {}
    for path in sorted(base.rglob("*.nc")):
        rel = path.relative_to(base)
        if len(rel.parts) < 4:
            continue
        gcm, rcm, feature_id = rel.parts[0:3]
        if allowed_gcms is not None and gcm not in allowed_gcms:
            continue
        if requested_aliases is not None:
            rcm_aliases = {rcm, rcm.replace("NARCliM2-0-WRF412", "")}
            if not (rcm_aliases & requested_aliases):
                continue
        if allowed_features is not None and str(feature_id) not in allowed_features:
            continue

        key = MemberKey(gcm=gcm, rcm=rcm, feature_id=str(feature_id))
        result.setdefault(key, []).append(path)
    return result


def _open_daily_member(*, files: list[Path], variable: str) -> xr.DataArray:
    arrays = []
    for path in sorted(files):
        with xr.open_dataset(
            path,
            engine="netcdf4",
            decode_times=True,
            decode_timedelta=False,
        ) as ds:
            var_name = _choose_variable(ds, variable)
            arrays.append(ds[var_name].load())

    if not arrays:
        raise RuntimeError("No arrays were opened.")

    combined = xr.concat(arrays, dim="time").sortby("time")
    time_index = combined.get_index("time")
    keep = ~time_index.duplicated()
    combined = combined.isel(time=np.flatnonzero(keep))

    if variable == "prAdjust":
        combined = _precipitation_to_mm_per_day(combined)
    else:
        combined = _temperature_to_celsius(combined, variable=variable)

    combined.name = variable
    return combined


def _derive_tasmean(tasmax: xr.DataArray, tasmin: xr.DataArray) -> xr.DataArray:
    tmax, tmin = xr.align(tasmax, tasmin, join="inner", copy=False)
    if tmax.sizes.get("time", 0) == 0:
        raise ValueError("tasmaxAdjust and tasminAdjust have no overlapping time steps.")
    out = (tmax + tmin) / 2.0
    out.name = "tasmeanAdjust"
    out.attrs = {
        "long_name": "Daily mean air temperature derived from daily max and min",
        "description": "tasmeanAdjust = (tasmaxAdjust + tasminAdjust) / 2",
        "units": "degC",
        "derived_from": "tasmaxAdjust,tasminAdjust",
    }
    return out


def _monthly_aggregate(da: xr.DataArray, *, variable: str) -> xr.DataArray:
    r = da.resample(time="MS")
    if variable == "prAdjust":
        out = r.sum(dim="time", skipna=True, min_count=1, keep_attrs=True)
        out.attrs["aggregation"] = "monthly total precipitation"
        out.attrs["units"] = "mm"
    else:
        out = r.mean(dim="time", skipna=True, keep_attrs=True)
        out.attrs["aggregation"] = "monthly mean temperature"
        out.attrs["units"] = "degC"
    out.name = variable
    return out


def _annual_aggregate(da: xr.DataArray, *, variable: str) -> xr.DataArray:
    r = da.resample(time="YS")
    if variable == "prAdjust":
        out = r.sum(dim="time", skipna=True, min_count=1, keep_attrs=True)
        out.attrs["aggregation"] = "annual total precipitation"
        out.attrs["units"] = "mm"
    else:
        out = r.mean(dim="time", skipna=True, keep_attrs=True)
        out.attrs["aggregation"] = "annual mean temperature"
        out.attrs["units"] = "degC"

    years = np.asarray(da["time"].dt.year.values, dtype=int)
    months = np.asarray(da["time"].dt.month.values, dtype=int)
    complete_years = []
    for year in np.unique(years):
        if set(months[years == year].tolist()) == set(range(1, 13)):
            complete_years.append(int(year))

    out_years = np.asarray(out["time"].dt.year.values, dtype=int)
    out = out.isel(time=np.flatnonzero(np.isin(out_years, complete_years)))
    out.name = variable
    return out


def _seasonal_aggregate(da: xr.DataArray, *, variable: str) -> xr.DataArray:
    month = np.asarray(da["time"].dt.month.values, dtype=int)
    year = np.asarray(da["time"].dt.year.values, dtype=int)

    season = np.empty(month.shape, dtype=object)
    season[np.isin(month, [12, 1, 2])] = "DJF"
    season[np.isin(month, [3, 4, 5])] = "MAM"
    season[np.isin(month, [6, 7, 8])] = "JJA"
    season[np.isin(month, [9, 10, 11])] = "SON"

    season_year = year.copy()
    season_year[month == 12] += 1
    group_keys = np.array(
        [sy * 10 + SEASON_ORDER[s] for sy, s in zip(season_year, season)],
        dtype=int,
    )

    work = da.assign_coords(
        _group=("time", group_keys),
        _season=("time", season.astype(str)),
        _season_year=("time", season_year),
        _month=("time", month),
    )

    expected = {
        "DJF": {12, 1, 2},
        "MAM": {3, 4, 5},
        "JJA": {6, 7, 8},
        "SON": {9, 10, 11},
    }

    outputs, times, seasons, years_out = [], [], [], []
    for key in np.unique(group_keys):
        part = work.where(work["_group"] == key, drop=True)
        season_name = str(part["_season"].values[0])
        sy = int(part["_season_year"].values[0])
        months_present = set(np.asarray(part["_month"].values, dtype=int).tolist())
        if months_present != expected[season_name]:
            continue

        if variable == "prAdjust":
            value = part.sum(dim="time", skipna=True, min_count=1, keep_attrs=True)
        else:
            value = part.mean(dim="time", skipna=True, keep_attrs=True)

        for name in ("_group", "_season", "_season_year", "_month"):
            if name in value.coords:
                value = value.drop_vars(name)

        outputs.append(value)
        seasons.append(season_name)
        years_out.append(sy)
        times.append(np.datetime64(
            f"{sy:04d}-{SEASON_REPRESENTATIVE_MONTH[season_name]:02d}-15"
        ))

    if not outputs:
        return da.isel(time=slice(0, 0))

    out = xr.concat(outputs, dim=pd.Index(times, name="time"))
    out = out.assign_coords(
        season=("time", np.asarray(seasons, dtype=str)),
        season_year=("time", np.asarray(years_out, dtype=int)),
    )
    if variable == "prAdjust":
        out.attrs["aggregation"] = "seasonal total precipitation"
        out.attrs["units"] = "mm"
    else:
        out.attrs["aggregation"] = "seasonal mean temperature"
        out.attrs["units"] = "degC"
    out.name = variable
    return out


def _aggregate(da: xr.DataArray, *, variable: str, frequency: str) -> xr.DataArray:
    if frequency == "Monthly":
        return _monthly_aggregate(da, variable=variable)
    if frequency == "Seasonal":
        return _seasonal_aggregate(da, variable=variable)
    if frequency == "Annual":
        return _annual_aggregate(da, variable=variable)
    raise ValueError(f"Unknown frequency {frequency!r}")


def _build_ensemble(
    member_series: dict[str, xr.DataArray],
    *,
    frequency: str,
) -> xr.DataArray:
    """
    Combine temporally aggregated ensemble-member time series.

    NARCliM ensemble members may use different source calendar/date objects.
    After aggregation, the scientifically relevant key is the analysis period
    (year-month, season-year/season, or year), not the original model-specific
    timestamp.

    To make the ensemble alignment robust AND keep NetCDF/CSV/plot time values
    standards-friendly, each member is assigned a synthetic Gregorian
    ``datetime64[ns]`` coordinate representing the same analysis period:

        Monthly  -> first day of month (YYYY-MM-01)
        Seasonal -> representative mid-season date
                    DJF=Feb-15, MAM=May-15, JJA=Aug-15, SON=Nov-15
        Annual   -> July 1 of the year

    The original seasonal meaning is also retained in ``season`` and
    ``season_year`` coordinates.
    """

    if not member_series:
        raise ValueError("No ensemble member series were supplied.")

    labels: list[str] = []
    arrays: list[xr.DataArray] = []

    for member_label, da in member_series.items():
        if "time" not in da.dims:
            raise ValueError(f"{member_label} has no time dimension.")

        # -------------------------------------------------------------
        # Build a common Gregorian analysis-time axis.
        # -------------------------------------------------------------
        if frequency == "Monthly":
            years = np.asarray(da["time"].dt.year.values, dtype=int)
            months = np.asarray(da["time"].dt.month.values, dtype=int)

            analysis_time = np.asarray(
                [np.datetime64(f"{year:04d}-{month:02d}-01")
                 for year, month in zip(years, months)],
                dtype="datetime64[ns]",
            )

            period_labels = np.asarray(
                [f"{year:04d}-{month:02d}"
                 for year, month in zip(years, months)],
                dtype=str,
            )

        elif frequency == "Seasonal":
            if "season_year" in da.coords and "season" in da.coords:
                years = np.asarray(da["season_year"].values, dtype=int)
                seasons = np.asarray(da["season"].values, dtype=str)
            else:
                years = np.asarray(da["time"].dt.year.values, dtype=int)
                months = np.asarray(da["time"].dt.month.values, dtype=int)
                month_to_season = {2: "DJF", 5: "MAM", 8: "JJA", 11: "SON"}
                try:
                    seasons = np.asarray(
                        [month_to_season[int(month)] for month in months],
                        dtype=str,
                    )
                except KeyError as exc:
                    raise ValueError(
                        "Seasonal series contains a representative month other "
                        f"than 2, 5, 8 or 11: {exc}"
                    ) from exc

            invalid_seasons = sorted(set(seasons) - set(SEASON_REPRESENTATIVE_MONTH))
            if invalid_seasons:
                raise ValueError(
                    f"{member_label}: unsupported seasons {invalid_seasons}."
                )

            analysis_time = np.asarray(
                [
                    np.datetime64(
                        f"{year:04d}-{SEASON_REPRESENTATIVE_MONTH[season]:02d}-15"
                    )
                    for year, season in zip(years, seasons)
                ],
                dtype="datetime64[ns]",
            )

            period_labels = np.asarray(
                [f"{year:04d}-{season}" for year, season in zip(years, seasons)],
                dtype=str,
            )

            # Reassign them explicitly so every member carries identical
            # coordinate types before alignment/concatenation.
            da = da.assign_coords(
                season=("time", seasons.astype(str)),
                season_year=("time", years.astype(int)),
            )

        elif frequency == "Annual":
            years = np.asarray(da["time"].dt.year.values, dtype=int)

            analysis_time = np.asarray(
                [np.datetime64(f"{year:04d}-07-01") for year in years],
                dtype="datetime64[ns]",
            )

            period_labels = np.asarray(
                [f"{year:04d}" for year in years],
                dtype=str,
            )

        else:
            raise ValueError(f"Unsupported frequency: {frequency}")

        da = da.assign_coords(
            time=("time", analysis_time),
            analysis_period=("time", period_labels),
        )

        # One value per analysis period per member is required.
        period_index = pd.Index(period_labels)
        if period_index.has_duplicates:
            duplicates = period_index[period_index.duplicated()].unique().tolist()
            raise ValueError(
                f"{member_label}: duplicate {frequency} periods found: "
                f"{duplicates[:10]}"
            )

        labels.append(member_label)
        arrays.append(da)

    # -------------------------------------------------------------
    # Align by the synthetic time axis and spatial coordinates.
    # -------------------------------------------------------------
    aligned = xr.align(*arrays, join="inner", copy=False)

    if not aligned or aligned[0].sizes.get("time", 0) == 0:
        periods_by_member = {
            label: (
                str(array["analysis_period"].values[0])
                if array.sizes.get("time", 0) else "EMPTY",
                str(array["analysis_period"].values[-1])
                if array.sizes.get("time", 0) else "EMPTY",
                int(array.sizes.get("time", 0)),
            )
            for label, array in zip(labels, arrays)
        }
        raise ValueError(
            "No common analysis periods remain across ensemble members.\n"
            f"Member period ranges: {periods_by_member}"
        )

    ensemble = xr.concat(
        aligned,
        dim=pd.Index(labels, name="ensemble_member"),
    )

    return ensemble

def _ensemble_statistics(ensemble: xr.DataArray) -> dict[str, xr.DataArray]:
    return {
        "min": ensemble.min("ensemble_member", skipna=True, keep_attrs=True),
        "mean": ensemble.mean("ensemble_member", skipna=True, keep_attrs=True),
        "max": ensemble.max("ensemble_member", skipna=True, keep_attrs=True),
    }


def _study_area_member_series(ensemble: xr.DataArray) -> xr.DataArray:
    spatial = _spatial_dims(ensemble)
    if not spatial:
        return ensemble
    return ensemble.mean(dim=spatial, skipna=True, keep_attrs=True)


def _area_ensemble_statistics(ensemble: xr.DataArray) -> dict[str, xr.DataArray]:
    member_area = _study_area_member_series(ensemble)
    return {
        "min": member_area.min("ensemble_member", skipna=True, keep_attrs=True),
        "mean": member_area.mean("ensemble_member", skipna=True, keep_attrs=True),
        "max": member_area.max("ensemble_member", skipna=True, keep_attrs=True),
    }


def _save_ensemble_netcdf(
    *, stat_da: xr.DataArray, variable: str, statistic: str, frequency: str,
    scenario: str, window_name: str, start_year: int, end_year: int,
    feature_id: str, n_members: int, output_path: Path, overwrite: bool,
) -> str:
    if output_path.exists() and not overwrite:
        return "skipped_existing"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    name = f"{variable}_{statistic}"
    out = stat_da.astype("float32")
    out.name = name
    out.attrs = stat_da.attrs.copy()
    out.attrs.update({
        "source_variable": variable,
        "ensemble_statistic": statistic,
        "temporal_frequency": frequency.lower(),
        "scenario": scenario,
        "time_horizon": window_name,
        "requested_start_year": int(start_year),
        "requested_end_year": int(end_year),
        "feature_id": str(feature_id),
        "ensemble_member_count": int(n_members),
        "analysis_type": "time-series ensemble uncertainty; no full-horizon temporal averaging",
    })
    ds = out.to_dataset()

    # Add useful metadata to the common analysis coordinates.
    if "analysis_period" in ds.coords:
        ds["analysis_period"].attrs.update({
            "long_name": "analysis period label",
            "description": (
                "Common period key used to align climate-model ensemble members "
                "independently of their original calendar encoding."
            ),
        })
    if "season" in ds.coords:
        ds["season"].attrs["long_name"] = "climatological season"
    if "season_year" in ds.coords:
        ds["season_year"].attrs["long_name"] = "year assigned to climatological season"
    ds.attrs.update({
        "title": f"NARCliM {frequency.lower()} {variable} ensemble {statistic} time series",
        "scenario": scenario,
        "time_horizon": window_name,
        "start_year": int(start_year),
        "end_year": int(end_year),
        "feature_id": str(feature_id),
        "ensemble_members": int(n_members),
    })
    temp = output_path.with_suffix(".tmp.nc")
    temp.unlink(missing_ok=True)
    ds.to_netcdf(
        temp, engine="netcdf4", mode="w",
        encoding={name: {"zlib": True, "complevel": 4, "dtype": "float32"}},
    )
    temp.replace(output_path)
    return "saved"


def _build_csv(
    *, area_stats: dict[str, xr.DataArray], variable: str, frequency: str,
    scenario: str, window_name: str, start_year: int, end_year: int,
    feature_id: str, n_members: int,
) -> pd.DataFrame:
    """Build the study-area ensemble time-series CSV.

    ``_build_ensemble`` guarantees a Gregorian datetime64 time coordinate, so
    CSV and plotting code never need to parse strings such as ``1985-MAM``.
    """
    mean_da = area_stats["mean"]

    time_values = pd.DatetimeIndex(mean_da["time"].values)

    df = pd.DataFrame({
        "time": time_values,
        "min": np.asarray(area_stats["min"].values, dtype=float),
        "mean": np.asarray(area_stats["mean"].values, dtype=float),
        "max": np.asarray(area_stats["max"].values, dtype=float),
    })

    if "analysis_period" in mean_da.coords:
        df.insert(1, "analysis_period", mean_da["analysis_period"].values.astype(str))

    if frequency == "Monthly":
        insert_at = 2 if "analysis_period" in df.columns else 1
        df.insert(insert_at, "year", df["time"].dt.year)
        df.insert(insert_at + 1, "month", df["time"].dt.month)
        df.insert(insert_at + 2, "month_name", df["time"].dt.month_name().str[:3])

    elif frequency == "Seasonal":
        if "season" not in mean_da.coords or "season_year" not in mean_da.coords:
            raise ValueError(
                "Seasonal ensemble output is missing season/season_year coordinates."
            )

        insert_at = 2 if "analysis_period" in df.columns else 1
        df.insert(
            insert_at,
            "season_year",
            np.asarray(mean_da["season_year"].values, dtype=int),
        )
        df.insert(
            insert_at + 1,
            "season",
            np.asarray(mean_da["season"].values, dtype=str),
        )

    elif frequency == "Annual":
        insert_at = 2 if "analysis_period" in df.columns else 1
        df.insert(insert_at, "year", df["time"].dt.year)

    df["variable"] = variable
    df["units"] = VARIABLE_UNITS[variable]
    df["scenario"] = scenario
    df["time_horizon"] = window_name
    df["start_year"] = int(start_year)
    df["end_year"] = int(end_year)
    df["feature_id"] = str(feature_id)
    df["ensemble_members"] = int(n_members)
    return df

def _plot_uncertainty(
    *, df: pd.DataFrame, member_area: xr.DataArray,
    variable: str, frequency: str, scenario: str,
    window_name: str, start_year: int, end_year: int, feature_id: str,
    output_path: Path, dpi: int, show: bool, overwrite: bool,
) -> str:
    """
    Plot all study-area ensemble-member time series together with the
    ensemble min-max envelope and ensemble mean.

    Individual GCM x RCM member series are drawn as thin light-grey lines
    and are intentionally excluded from the legend.
    """
    if output_path.exists() and not overwrite:
        return "skipped_existing"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 5.5))

    x = pd.DatetimeIndex(df["time"])
    y_min = df["min"].to_numpy(dtype=float)
    y_mean = df["mean"].to_numpy(dtype=float)
    y_max = df["max"].to_numpy(dtype=float)

    # ==============================================================
    # INDIVIDUAL ENSEMBLE MEMBERS
    #
    # Plot each GCM × RCM member with a different color.
    # Model names are intentionally excluded from the legend.
    # ==============================================================

    member_times = pd.DatetimeIndex(
        member_area["time"].values
    )

    members = member_area["ensemble_member"].values

    # tab10 provides 10 clearly distinguishable colors,
    # which matches the usual 5 GCM × 2 RCM = 10 members.
    cmap = plt.get_cmap("tab10")

    for i, member in enumerate(members):

        member_values = np.asarray(
            member_area.sel(
                ensemble_member=member
            ).values,
            dtype=float,
        )

        ax.plot(
            member_times,
            member_values,
            color=cmap(i % 10),
            linewidth=0.8,
            alpha=0.65,
            label="_nolegend_",
            zorder=1,
        )

    # ==============================================================
    # ENSEMBLE MIN-MAX RANGE
    # ==============================================================

    ax.fill_between(
        x,
        y_min,
        y_max,
        color="steelblue",
        alpha=0.25,
        label="Min-max range",
        zorder=2,
    )

    # ==============================================================
    # ENSEMBLE MEAN
    # ==============================================================

    ax.plot(
        x,
        y_mean,
        color="navy",
        linewidth=1.8,
        label="Mean",
        zorder=4,
    )

    # ==============================================================
    # MINIMUM / MAXIMUM ENVELOPE BOUNDARIES
    # ============================================================== 

    ax.plot(
        x,
        y_min,
        color="steelblue",
        linewidth=0.7,
        alpha=0.75,
        label="_nolegend_",
        zorder=3,
    )

    ax.plot(
        x,
        y_max,
        color="steelblue",
        linewidth=0.7,
        alpha=0.75,
        label="_nolegend_",
        zorder=3,
    )

    label = VARIABLE_LABELS[variable]
    units = VARIABLE_UNITS[variable]

    temporal_text = (
        {
            "Monthly": "monthly total",
            "Seasonal": "seasonal total",
            "Annual": "annual total",
        }[frequency]
        if variable == "prAdjust"
        else
        {
            "Monthly": "monthly mean",
            "Seasonal": "seasonal mean",
            "Annual": "annual mean",
        }[frequency]
    )

    ax.set_title(
        f"{label} – {temporal_text} - {scenario} ({start_year}-{end_year})"
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(f"{label} ({units})")
    ax.grid(True, alpha=0.25)

    # Only the summary elements appear in the legend.
    ax.legend(frameon=False, ncol=1)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return "saved"


def _load_variable_members(
    *, data_root: Path, variable: str, scenario: str, window_name: str,
    gcms: Iterable[str] | None, rcms: Iterable[str] | None,
    feature_ids: Iterable[str] | None, verbose: bool,
) -> dict[MemberKey, xr.DataArray]:
    discovered = _discover_member_files(
        data_root=data_root, variable=variable, scenario=scenario, window=window_name,
        gcms=gcms, rcms=rcms, feature_ids=feature_ids,
    )
    out = {}
    for key, paths in discovered.items():
        _log(
            f"    [OPEN] {variable} | {key.gcm} | {key.rcm} | "
            f"{key.feature_id} | files={len(paths)}",
            verbose=verbose,
        )
        out[key] = _open_daily_member(files=paths, variable=variable)
    return out


def _prepare_daily_variables(
    *, data_root: Path, scenario: str, window_name: str,
    gcms: Iterable[str] | None, rcms: Iterable[str] | None,
    feature_ids: Iterable[str] | None, variables: Iterable[str], verbose: bool,
) -> dict[str, dict[MemberKey, xr.DataArray]]:
    requested = set(variables)
    invalid = requested - set(ALL_VARIABLES)
    if invalid:
        raise ValueError(f"Unsupported variables: {sorted(invalid)}")

    needed_sources = set(requested) & set(SOURCE_VARIABLES)
    if "tasmeanAdjust" in requested:
        needed_sources.update({"tasmaxAdjust", "tasminAdjust"})

    daily: dict[str, dict[MemberKey, xr.DataArray]] = {}
    for variable in SOURCE_VARIABLES:
        if variable in needed_sources:
            daily[variable] = _load_variable_members(
                data_root=data_root, variable=variable, scenario=scenario,
                window_name=window_name, gcms=gcms, rcms=rcms,
                feature_ids=feature_ids, verbose=verbose,
            )

    if "tasmeanAdjust" in requested:
        tmax = daily.get("tasmaxAdjust", {})
        tmin = daily.get("tasminAdjust", {})
        common = set(tmax) & set(tmin)
        if set(tmax) - set(tmin):
            warnings.warn("Some tasmaxAdjust members have no matching tasminAdjust member.")
        if set(tmin) - set(tmax):
            warnings.warn("Some tasminAdjust members have no matching tasmaxAdjust member.")
        daily["tasmeanAdjust"] = {
            key: _derive_tasmean(tmax[key], tmin[key]) for key in common
        }

    return {v: daily[v] for v in requested if v in daily}

def _discover_scenarios_and_time_windows(
    *,
    data_root: Path,
    variables: Iterable[str],
) -> tuple[tuple[str, ...], dict[str, tuple[int, int]]]:
    """
    Discover scenarios and time windows from data already downloaded
    by workflow.py.

    Expected structure
    ------------------
    NARCliM Data/
        variable/
            scenario/
                time_window/
                    GCM/
                        RCM/
                            feature_id/
                                *.nc

    The actual start/end years are read from the NetCDF metadata written
    by workflow.py. If those attributes are unavailable, the years are
    determined from the NetCDF time coordinate.
    """

    scenarios_found = set()
    windows_found = {}

    # -----------------------------------------------------------------
    # tasmeanAdjust is derived locally and therefore does not have its
    # own downloaded source folder.
    # -----------------------------------------------------------------

    source_variables = []

    for variable in variables:

        if variable == "tasmeanAdjust":

            for source in (
                "tasmaxAdjust",
                "tasminAdjust",
            ):
                if source not in source_variables:
                    source_variables.append(source)

        else:

            if variable not in source_variables:
                source_variables.append(variable)

    # -----------------------------------------------------------------
    # Search downloaded variable folders.
    # -----------------------------------------------------------------

    for variable in source_variables:

        variable_root = (
            data_root
            / variable
        )

        if not variable_root.exists():
            continue

        # -------------------------------------------------------------
        # Scenario folders
        # -------------------------------------------------------------

        for scenario_dir in sorted(
            path
            for path in variable_root.iterdir()
            if path.is_dir()
        ):

            scenario = scenario_dir.name

            scenarios_found.add(
                scenario
            )

            # ---------------------------------------------------------
            # Time-window folders
            # ---------------------------------------------------------

            for window_dir in sorted(
                path
                for path in scenario_dir.iterdir()
                if path.is_dir()
            ):

                window_name = (
                    window_dir.name
                )

                # Find one downloaded NetCDF belonging to this window.
                example_file = next(
                    window_dir.rglob("*.nc"),
                    None,
                )

                if example_file is None:
                    continue

                # -----------------------------------------------------
                # Read the requested years from workflow.py metadata.
                # -----------------------------------------------------

                with xr.open_dataset(
                    example_file,
                    engine="netcdf4",
                    decode_times=True,
                    decode_timedelta=False,
                ) as ds:

                    start_year = ds.attrs.get(
                        "requested_start_year"
                    )

                    end_year = ds.attrs.get(
                        "requested_end_year"
                    )

                    # -------------------------------------------------
                    # Fallback:
                    # derive years from actual NetCDF time coordinate.
                    # -------------------------------------------------

                    if (
                        start_year is None
                        or end_year is None
                    ):

                        if "time" not in ds.coords:

                            raise ValueError(
                                "Could not determine the time "
                                f"window for:\n{example_file}\n"
                                "The file contains neither "
                                "requested_start_year / "
                                "requested_end_year metadata "
                                "nor a time coordinate."
                            )

                        years = np.asarray(
                            ds["time"]
                            .dt.year
                            .values,
                            dtype=int,
                        )

                        if years.size == 0:

                            raise ValueError(
                                "Empty time coordinate in:\n"
                                f"{example_file}"
                            )

                        start_year = int(
                            years.min()
                        )

                        end_year = int(
                            years.max()
                        )

                year_range = (
                    int(start_year),
                    int(end_year),
                )

                # -----------------------------------------------------
                # Check that a named time window has a consistent
                # definition across variables/scenarios.
                # -----------------------------------------------------

                if window_name in windows_found:

                    existing_range = (
                        windows_found[
                            window_name
                        ]
                    )

                    if existing_range != year_range:

                        raise ValueError(
                            "Conflicting year ranges were "
                            f"found for '{window_name}'.\n"
                            f"Existing: {existing_range}\n"
                            f"Found:    {year_range}\n"
                            f"File: {example_file}"
                        )

                else:

                    windows_found[
                        window_name
                    ] = year_range

    # -----------------------------------------------------------------
    # Validation
    # -----------------------------------------------------------------

    if not scenarios_found:

        raise FileNotFoundError(
            "No downloaded NARCliM scenarios "
            "were found under:\n"
            f"{data_root}"
        )

    if not windows_found:

        raise FileNotFoundError(
            "No downloaded NARCliM time windows "
            "containing NetCDF files were found under:\n"
            f"{data_root}"
        )

    # Historical first, followed by SSP scenarios.
    scenarios_found = tuple(
        sorted(
            scenarios_found,
            key=lambda value: (
                0
                if value == "historical"
                else 1,
                value,
            ),
        )
    )

    return (
        scenarios_found,
        windows_found,
    )

def run_rain_temp_timeseries(
    *, input_root: str | Path,
    variables: Iterable[str] = ALL_VARIABLES,
    frequencies: Iterable[str] = FREQUENCIES,
    scenarios: Iterable[str] | None = None,
    time_windows: dict[str, tuple[int, int]] | None = None,
    gcms: Iterable[str] | None = None,
    rcms: Iterable[str] | None = None,
    feature_ids: Iterable[str] | None = None,
    output_folder_name: str = "Climate_Indicies",
    overwrite: bool = False,
    plot_dpi: int = 300,
    show_plots: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Calculate monthly, seasonal and annual ensemble time series for
    prAdjust, tasmaxAdjust, tasminAdjust and derived tasmeanAdjust.

    By default, scenarios and time windows are discovered automatically from
    the NetCDF subsets already downloaded by workflow.py. They can still be
    supplied explicitly as optional filters/overrides.

    Rainfall uses temporal totals; temperature uses temporal means.
    The time dimension is preserved. This is not a whole-horizon temporal mean.
    """
    project_root, data_root = _normalise_root(input_root)
    variables = tuple(variables)
    frequencies = tuple(f.title() for f in frequencies)

    unknown = set(frequencies) - set(FREQUENCIES)
    if unknown:
        raise ValueError(f"Unknown frequencies: {sorted(unknown)}")

    discovered_scenarios, discovered_windows = (
        _discover_scenarios_and_time_windows(
            data_root=data_root,
            variables=variables,
        )
    )

    if scenarios is None:
        scenarios = discovered_scenarios
    else:
        scenarios = tuple(scenarios)

    if time_windows is None:
        time_windows = discovered_windows
    else:
        time_windows = dict(time_windows)

    output_root = project_root / output_folder_name / "Rain_Temp_Timeseries"
    output_root.mkdir(parents=True, exist_ok=True)
    records = []

    _log("=" * 100, verbose=verbose)
    _log("NARCliM RAINFALL / TEMPERATURE TIME-SERIES WORKFLOW", verbose=verbose)
    _log("=" * 100, verbose=verbose)
    _log(f"[DISCOVERED SCENARIOS] {discovered_scenarios}", verbose=verbose)
    _log(f"[DISCOVERED WINDOWS] {discovered_windows}", verbose=verbose)
    _log(f"[PROCESSING SCENARIOS] {tuple(scenarios)}", verbose=verbose)
    _log(f"[PROCESSING WINDOWS] {time_windows}", verbose=verbose)


    for scenario in scenarios:
        for window_name, (start_year, end_year) in time_windows.items():
            if scenario == "historical" and start_year > 2014:
                continue
            if scenario != "historical" and end_year < 2015:
                continue

            _log(
                f"\n[GROUP] {scenario} | {window_name} | {start_year}-{end_year}",
                verbose=verbose,
            )

            daily_by_variable = _prepare_daily_variables(
                data_root=data_root, scenario=scenario, window_name=window_name,
                gcms=gcms, rcms=rcms, feature_ids=feature_ids,
                variables=variables, verbose=verbose,
            )

            for variable in variables:
                member_daily = daily_by_variable.get(variable, {})
                if not member_daily:
                    records.append({
                        "status": "no_data", "variable": variable,
                        "scenario": scenario, "time_window": window_name,
                        "frequency": None, "feature_id": None,
                        "output_type": None, "output": None,
                        "message": "No matching local members found.",
                    })
                    continue

                available_features = sorted({key.feature_id for key in member_daily})
                for feature_id in available_features:
                    feature_members = {
                        key.label: da for key, da in member_daily.items()
                        if key.feature_id == feature_id
                    }
                    for frequency in frequencies:
                        aggregated_members = {}
                        for label, daily_da in feature_members.items():
                            agg = _aggregate(daily_da, variable=variable, frequency=frequency)
                            if agg.sizes.get("time", 0) > 0:
                                aggregated_members[label] = agg

                        if not aggregated_members:
                            continue

                        ensemble = _build_ensemble(
                            aggregated_members,
                            frequency=frequency,
                        )

                        grid_stats = _ensemble_statistics(
                            ensemble
                        )

                        # Spatially averaged series for each individual
                        # GCM x RCM member. These are used for the light-grey
                        # model lines on the plot.
                        member_area = _study_area_member_series(
                            ensemble
                        )

                        area_stats = _area_ensemble_statistics(
                            ensemble
                        )

                        base = output_root / scenario / window_name / frequency

                        for statistic in STATISTICS:
                            nc_name = (
                                f"{frequency}_{variable}_{statistic}_{scenario}_"
                                f"{start_year}_{end_year}_{_safe_name(feature_id)}.nc"
                            )
                            nc_path = base / "NetCDF" / variable / nc_name
                            status = _save_ensemble_netcdf(
                                stat_da=grid_stats[statistic], variable=variable,
                                statistic=statistic, frequency=frequency,
                                scenario=scenario, window_name=window_name,
                                start_year=start_year, end_year=end_year,
                                feature_id=feature_id,
                                n_members=len(aggregated_members),
                                output_path=nc_path, overwrite=overwrite,
                            )
                            records.append({
                                "status": status, "variable": variable,
                                "scenario": scenario, "time_window": window_name,
                                "frequency": frequency, "feature_id": feature_id,
                                "ensemble_statistic": statistic,
                                "ensemble_members": len(aggregated_members),
                                "output_type": "NetCDF", "output": str(nc_path),
                                "message": "",
                            })

                        friendly = {
                            "prAdjust": "rainfall",
                            "tasmaxAdjust": "tasmax",
                            "tasminAdjust": "tasmin",
                            "tasmeanAdjust": "tasmean",
                        }[variable]
                        csv_path = (
                            base / "CSV" / variable /
                            f"{frequency}_{friendly}_{scenario}_{start_year}_{end_year}_"
                            f"{_safe_name(feature_id)}.csv"
                        )
                        csv_path.parent.mkdir(parents=True, exist_ok=True)
                        df = _build_csv(
                            area_stats=area_stats, variable=variable,
                            frequency=frequency, scenario=scenario,
                            window_name=window_name, start_year=start_year,
                            end_year=end_year, feature_id=feature_id,
                            n_members=len(aggregated_members),
                        )
                        if csv_path.exists() and not overwrite:
                            csv_status = "skipped_existing"
                        else:
                            df.to_csv(csv_path, index=False)
                            csv_status = "saved"
                        records.append({
                            "status": csv_status, "variable": variable,
                            "scenario": scenario, "time_window": window_name,
                            "frequency": frequency, "feature_id": feature_id,
                            "ensemble_statistic": "min_mean_max",
                            "ensemble_members": len(aggregated_members),
                            "output_type": "CSV", "output": str(csv_path),
                            "message": "",
                        })

                        plot_friendly = {
                            "prAdjust": "rainfall",
                            "tasmaxAdjust": "maximum_temperature",
                            "tasminAdjust": "minimum_temperature",
                            "tasmeanAdjust": "mean_temperature",
                        }[variable]
                        png_path = (
                            base / "Plots" / variable /
                            f"{frequency}_{plot_friendly}_uncertainty_{scenario}_"
                            f"{start_year}_{end_year}_{_safe_name(feature_id)}.png"
                        )
                        png_status = _plot_uncertainty(
                            df=df,
                            member_area=member_area,
                            variable=variable,
                            frequency=frequency,
                            scenario=scenario,
                            window_name=window_name,
                            start_year=start_year,
                            end_year=end_year,
                            feature_id=feature_id,
                            output_path=png_path,
                            dpi=plot_dpi,
                            show=show_plots,
                            overwrite=overwrite,
                        )
                        records.append({
                            "status": png_status, "variable": variable,
                            "scenario": scenario, "time_window": window_name,
                            "frequency": frequency, "feature_id": feature_id,
                            "ensemble_statistic": "min_mean_max",
                            "ensemble_members": len(aggregated_members),
                            "output_type": "PNG", "output": str(png_path),
                            "message": "",
                        })

                        _log(
                            f"[DONE] {variable} | {scenario} | {window_name} | "
                            f"{feature_id} | {frequency} | members={len(aggregated_members)}",
                            verbose=verbose,
                        )

    manifest = pd.DataFrame(records)
    manifest_path = output_root / "Rain_Temp_Timeseries_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    _log(f"\n[MANIFEST] {manifest_path}", verbose=verbose)
    return manifest