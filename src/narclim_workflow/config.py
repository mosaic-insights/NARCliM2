from __future__ import annotations

from dataclasses import dataclass


# =============================================================================
# THREDDS CONFIGURATION
# =============================================================================

THREDDS_ROOT = "https://thredds.nci.org.au"

PROJECT_ID = "zz63"

# -------------------------------------------------------------------------
# Default NARCliM collection.
#
# Most climate indices and bias-adjusted variables used by this package
# come from:
#
#   NARCliM2-0-derived/output-CMIP6/
#
# Some raw NARCliM variables, such as sfcWindmax, instead come from:
#
#   NARCliM2-0/output-CMIP6/
#
# The VariableSpec.collection field below allows this to be configured
# separately for each variable.
# -------------------------------------------------------------------------

COLLECTION = "NARCliM2-0-derived"

CMIP_OUTPUT = "output-CMIP6"

PROVIDER = "NSW-Government"


# =============================================================================
# DOMAINS, MODELS, AND SCENARIOS
# =============================================================================

DOMAINS = {
    "seaus_4km": "NARCliM2-0-SEAus-04",
    "aus_18": "AUS-18",
}


GCM_REALISATIONS = {
    "ACCESS-ESM1-5": "r6i1p1f1",
    "EC-Earth3-Veg": "r1i1p1f1",
    "MPI-ESM1-2-HR": "r1i1p1f1",
    "NorESM2-MM": "r1i1p1f1",
    "UKESM1-0-LL": "r1i1p1f2",
}


RCMS = (
    "NARCliM2-0-WRF412R3",
    "NARCliM2-0-WRF412R5",
)


SCENARIOS = (
    "historical",
    "ssp126",
    "ssp245",
    "ssp370",
)


# =============================================================================
# VARIABLE CONFIGURATION
# =============================================================================

@dataclass(frozen=True)
class VariableSpec:
    """
    Configuration describing where and how a NARCliM variable is stored.

    Parameters
    ----------
    name
        Variable name in the source NARCliM NetCDF file.

    output_type
        THREDDS product branch, for example:
            "bias-adjusted-output"
            "DD"

    processing_version
        NARCliM processing/version directory.

    frequency
        Source temporal frequency:
            "day"
            "mon"
            "yr"

    description
        Human-readable description shown by the workflow.

    typical_storage
        General source-file organisation, for example:
            "one_file_per_year"
            "whole_scenario_period"

    collection
        NCI NARCliM collection.

        Most variables use:
            NARCliM2-0-derived

        Raw model variables such as sfcWindmax use:
            NARCliM2-0

    available_scenarios
        Scenarios for which the variable is expected to be available.

    derived_variable
        Name of a derived package variable created from the source variable.

        None means the source variable is saved directly.

    threshold_kmh
        Optional wind-speed threshold used for derived storm indicators.
    """

    name: str
    output_type: str
    processing_version: str
    frequency: str
    description: str
    typical_storage: str

    # Most existing package variables are stored under NARCliM2-0-derived.
    collection: str = COLLECTION

    # Most variables are available for historical and all projections.
    available_scenarios: tuple[str, ...] = (
        "historical",
        "ssp126",
        "ssp245",
        "ssp370",
    )

    # ---------------------------------------------------------------------
    # Optional derived-variable configuration.
    #
    # These are mainly used for variables such as StormDaysGT89, where
    # workflow.py downloads one source variable (sfcWindmax) and derives
    # another variable from it.
    # ---------------------------------------------------------------------

    derived_variable: str | None = None

    threshold_kmh: float | None = None


# =============================================================================
# VARIABLES
# =============================================================================

VARIABLES = {

    # -------------------------------------------------------------------------
    # BIAS-ADJUSTED DAILY VARIABLES
    # -------------------------------------------------------------------------

    "prAdjust": VariableSpec(
        name="prAdjust",
        output_type="bias-adjusted-output",
        processing_version=(
            "v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009"
        ),
        frequency="day",
        description="Daily bias-adjusted precipitation.",
        typical_storage="one_file_per_year",
    ),


    "tasmaxAdjust": VariableSpec(
        name="tasmaxAdjust",
        output_type="bias-adjusted-output",
        processing_version=(
            "v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009"
        ),
        frequency="day",
        description="Daily bias-adjusted maximum temperature.",
        typical_storage="one_file_per_year",
    ),


    "tasminAdjust": VariableSpec(
        name="tasminAdjust",
        output_type="bias-adjusted-output",
        processing_version=(
            "v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009"
        ),
        frequency="day",
        description="Daily bias-adjusted minimum temperature.",
        typical_storage="one_file_per_year",
    ),


    # -------------------------------------------------------------------------
    # BIAS-ADJUSTED CLIMATE INDICES
    # -------------------------------------------------------------------------

    "TXge35": VariableSpec(
        name="TXge35",
        output_type="bias-adjusted-output",
        processing_version="v1-r1",
        frequency="yr",
        description=(
            "Yearly number of days with maximum temperature >=35 C."
        ),
        typical_storage="whole_scenario_period",
    ),


    "TNlt2": VariableSpec(
        name="TNlt2",
        output_type="bias-adjusted-output",
        processing_version="v1-r1",
        frequency="yr",
        description=(
            "Yearly number of days with minimum temperature <2 C."
        ),
        typical_storage="whole_scenario_period",
    ),


    # -------------------------------------------------------------------------
    # DERIVED CLIMATE INDICES
    # -------------------------------------------------------------------------

    "FFDIgt50": VariableSpec(
        name="FFDIgt50",
        output_type="DD",
        processing_version="v1-r1",
        frequency="yr",
        description=(
            "Yearly number of days with FFDI >=50."
        ),
        typical_storage="whole_scenario_period",
    ),


    "R20mm": VariableSpec(
        name="R20mm",
        output_type="DD",
        processing_version="v1-r1",
        frequency="yr",
        description=(
            "Yearly number of days with precipitation >=20 mm."
        ),
        typical_storage="whole_scenario_period",
    ),


    "R99p": VariableSpec(
        name="R99p",
        output_type="DD",
        processing_version="v1-r1",
        frequency="yr",
        description=(
            "Yearly total precipitation from extremely wet days."
        ),
        typical_storage="whole_scenario_period",
    ),


    "SPI12": VariableSpec(
        name="SPI12",
        output_type="DD",
        processing_version="v1-r1",
        frequency="mon",
        description=(
            "12-month Standardised Precipitation Index."
        ),
        typical_storage="whole_scenario_period",

        # SPI12 is not available in the historical experiment.
        available_scenarios=(
            "ssp126",
            "ssp245",
            "ssp370",
        ),
    ),


    # -------------------------------------------------------------------------
    # STORM HAZARD
    #
    # Source:
    #   NARCliM daily maximum near-surface wind speed (sfcWindmax)
    #
    # Source units:
    #   m s-1
    #
    # Derived package variable:
    #   StormDaysGT89
    #
    # Definition:
    #   Annual number of days where daily maximum near-surface wind speed
    #   exceeds 89 km/h.
    #
    #   89 km/h = 24.7222 m/s
    #
    # IMPORTANT:
    # sfcWindmax is stored under the RAW NARCliM collection:
    #
    #   NARCliM2-0/output-CMIP6/DD/
    #
    # rather than:
    #
    #   NARCliM2-0-derived/output-CMIP6/DD/
    # -------------------------------------------------------------------------

    "StormDays": VariableSpec(
        # Original NARCliM source variable.
        name="sfcWindmax",

        # Raw NARCliM DD collection.
        output_type="DD",

        processing_version="v1-r1",

        # Source data are daily.
        frequency="day",

        description=(
            "Annual number of days exceeding a user-defined "
            "daily maximum near-surface wind-speed threshold."
        ),

        typical_storage="one_file_per_year",

        # sfcWindmax comes from the raw NARCliM collection.
        collection="NARCliM2-0",

        available_scenarios=(
            "historical",
            "ssp126",
            "ssp245",
            "ssp370",
        ),

        # This only tells workflow.py that the source must be converted
        # to a storm-day exceedance index. The final variable name is
        # generated dynamically from the threshold.
        derived_variable="storm_days_threshold",

        # Default threshold.
        # The user can override this in run_workflow().
        threshold_kmh=89.0,
    ),
}