# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)

# flake8: noqa: F401
import importlib.metadata

from . import calibrations, conversions, corrections, normalize, reductions

try:
    __version__ = importlib.metadata.version(__package__ or __name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

providers = (
    *conversions.providers,
    *corrections.providers,
    *calibrations.providers,
    *normalize.providers,
    *reductions.providers,
)
"""
List of providers for setting up a Sciline pipeline.

This does not constitute a complete workflow but only
a skeleton for a generic reflectometry setting.
For an example of a complete workflow
see :py:data:`essreflectometry.amor.providers`.
"""

del importlib
