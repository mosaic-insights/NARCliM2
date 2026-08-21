from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
import re
import xml.etree.ElementTree as ET

import requests

from .config import (
    CMIP_OUTPUT,
    GCM_REALISATIONS,
    PROJECT_ID,
    PROVIDER,
    THREDDS_ROOT,
    VariableSpec,
)


# =============================================================================
# CATALOG DATA STRUCTURE
# =============================================================================

@dataclass(frozen=True)
class DatasetFile:
    """
    Representation of one NetCDF file discovered in a THREDDS catalog.
    """

    name: str
    url_path: str
    catalog_url: str

    @property
    def opendap_url(self) -> str:
        """
        Return the OPeNDAP URL used by xarray.
        """

        return (
            f"{THREDDS_ROOT}/thredds/dodsC/"
            f"{self.url_path}"
        )

    @property
    def file_server_url(self) -> str:
        """
        Return the direct THREDDS file-server URL.
        """

        return (
            f"{THREDDS_ROOT}/thredds/fileServer/"
            f"{self.url_path}"
        )


# =============================================================================
# BUILD THREDDS CATALOG URL
# =============================================================================

def build_catalog_url(
    *,
    domain: str,
    gcm: str,
    scenario: str,
    rcm: str,
    variable: VariableSpec,
    release: str = "latest",
) -> str:
    """
    Build the NCI THREDDS catalog URL for a NARCliM variable.

    The NARCliM collection is defined separately for each variable through
    VariableSpec.collection.

    Examples
    --------
    Most derived variables:

        NARCliM2-0-derived/output-CMIP6/...

    Raw model variables such as sfcWindmax:

        NARCliM2-0/output-CMIP6/...
    """

    realisation = GCM_REALISATIONS[
        gcm
    ]

    parts = [
        THREDDS_ROOT,
        "thredds",
        "catalog",
        PROJECT_ID,

        # IMPORTANT:
        # Use the collection configured for the individual variable.
        #
        # Examples:
        #   NARCliM2-0-derived
        #   NARCliM2-0
        variable.collection,

        CMIP_OUTPUT,
        variable.output_type,
        domain,
        PROVIDER,
        gcm,
        scenario,
        realisation,
        rcm,
        variable.processing_version,
        variable.frequency,
        variable.name,
        release,
        "catalog.xml",
    ]

    return "/".join(
        part.strip("/")
        for part in parts
    )


# =============================================================================
# READ THREDDS CATALOG
# =============================================================================

def read_catalog(
    catalog_url: str,
    timeout: int = 120,
) -> list[DatasetFile]:
    """
    Read a THREDDS catalog and return the NetCDF datasets it contains.
    """

    response = requests.get(
        catalog_url,
        timeout=timeout,
    )

    response.raise_for_status()

    root = ET.fromstring(
        response.content
    )

    ns = {
        "t": (
            "http://www.unidata.ucar.edu/"
            "namespaces/thredds/InvCatalog/v1.0"
        )
    }

    found = []

    for element in root.findall(
        ".//t:dataset",
        ns,
    ):

        url_path = (
            element.attrib.get(
                "urlPath"
            )
        )

        name = (
            element.attrib.get(
                "name",
                "",
            )
        )

        if (
            url_path
            and url_path.lower().endswith(
                (
                    ".nc",
                    ".nc4",
                )
            )
        ):

            found.append(
                DatasetFile(
                    name=(
                        name
                        or PurePosixPath(
                            url_path
                        ).name
                    ),
                    url_path=url_path,
                    catalog_url=catalog_url,
                )
            )

    if not found:
        raise RuntimeError(
            "No NetCDF files found in catalog: "
            f"{catalog_url}"
        )

    return sorted(
        found,
        key=lambda item: item.name,
    )


# =============================================================================
# INFER FILE YEAR RANGE
# =============================================================================

_DATE_RANGE_RE = re.compile(
    r"(?P<start>\d{4})"
    r"(?:\d{4})?"
    r"-"
    r"(?P<end>\d{4})"
    r"(?:\d{4})?"
    r"\.nc$",
    re.I,
)


def infer_file_year_range(
    filename: str,
) -> tuple[int, int] | None:
    """
    Infer the first and last year represented by a NetCDF filename.

    Supported examples
    ------------------
    Yearly daily files:

        sfcWindmax_..._19510101-19511231.nc
            -> (1951, 1951)

        prAdjust_..._19850101-19851231.nc
            -> (1985, 1985)

    Whole-period annual files:

        TXge35_..._1951-2014.nc
            -> (1951, 2014)

    Monthly ranges:

        SPI12_..._195101-210012.nc
            -> (1951, 2100)
    """

    match = _DATE_RANGE_RE.search(
        filename
    )

    if not match:
        return None

    return (
        int(
            match.group(
                "start"
            )
        ),
        int(
            match.group(
                "end"
            )
        ),
    )


# =============================================================================
# SELECT FILES OVERLAPPING A TIME WINDOW
# =============================================================================

def select_files_for_window(
    files: list[DatasetFile],
    start_year: int,
    end_year: int,
) -> list[DatasetFile]:
    """
    Return source NetCDF files whose represented years overlap the requested
    time window.

    Files whose date range cannot be inferred are retained as a fallback only
    when no date-parsable files were selected.
    """

    selected = []
    unknown = []

    for item in files:

        year_range = infer_file_year_range(
            item.name
        )

        if year_range is None:

            unknown.append(
                item
            )

        elif (
            year_range[0] <= end_year
            and year_range[1] >= start_year
        ):

            selected.append(
                item
            )

    return (
        selected
        or unknown
    )