"""
Public interface for the NARCliM2 workflow package.

The package provides workflows for:

- downloading and subsetting NARCliM2 climate data;
- ensemble uncertainty analysis;
- rainfall and temperature time-series analysis;
- approximate climate-adjusted flood hazard analysis.
"""

from .config import (
    GCM_REALISATIONS,
    RCMS,
    SCENARIOS,
    VARIABLES,
    VariableSpec,
)

from .workflow import run_workflow

from .uncertainty import (
    run_uncertainty_workflow,
)

from .rain_temp_timeseries import (
    run_rain_temp_timeseries,
)

from .flood import (
    run_flood_hazard,
)


__all__ = [
    # Configuration
    "GCM_REALISATIONS",
    "RCMS",
    "SCENARIOS",
    "VARIABLES",
    "VariableSpec",

    # Main workflows
    "run_workflow",
    "run_uncertainty_workflow",
    "run_rain_temp_timeseries",
    "run_flood_hazard",
]