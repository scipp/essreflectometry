# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
# flake8: noqa: F401
from itertools import chain

import scipp as sc

from . import beamline, calibrations, conversions, load, normalize, resolution, tools

# from .beamline import instrument_view_components
from .instrument_view import instrument_view
from .types import *

providers = list(
    chain(
        load.providers,
        calibrations.providers,
        conversions.providers,
        normalize.providers,
        resolution.providers,
        beamline.providers,
    )
)

default_parameters = {
    Supermirror[MValue]: sc.scalar(5, unit=sc.units.dimensionless),
    Supermirror[CriticalEdge]: 0.022 * sc.Unit('1/angstrom'),
    Supermirror[Alpha]: sc.scalar(0.25 / 0.088, unit=sc.units.angstrom),
}
