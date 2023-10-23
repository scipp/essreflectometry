# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
# flake8: noqa: F401
from itertools import chain

from . import calibrations, conversions, load, normalize, resolution, tools
from .beamline import instrument_view_components
from .instrument_view import instrument_view

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
