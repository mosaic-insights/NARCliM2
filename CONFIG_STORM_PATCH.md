# Required `config.py` storm entry

The configurable storm workflow expects a generic `StormDays` entry in
`VARIABLES`. Merge the following into your existing `config.py`.

Your `VariableSpec` dataclass should include these optional fields:

```python
collection: str = COLLECTION
derived_variable: str | None = None
threshold_kmh: float | None = None
```

Add:

```python
"StormDays": VariableSpec(
    name="sfcWindmax",
    output_type="DD",
    processing_version="v1-r1",
    frequency="day",
    description=(
        "Annual number of days exceeding a user-defined daily maximum "
        "near-surface wind-speed threshold."
    ),
    typical_storage="one_file_per_year",
    collection="NARCliM2-0",
    available_scenarios=(
        "historical",
        "ssp126",
        "ssp245",
        "ssp370",
    ),
    derived_variable="storm_days_threshold",
    threshold_kmh=89.0,
),
```

The `threshold_kmh=89.0` value is a default only. Users can override it with:

```python
run_workflow(
    ...,
    storm_threshold_kmh=63.0,
)
```

which creates `StormDaysGT63`.
