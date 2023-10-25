# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)

# flake8: noqa: F401
import importlib.metadata
import itertools

from . import conversions, corrections, io

try:
    __version__ = importlib.metadata.version(__package__ or __name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

providers = list(
    itertools.chain(
        conversions.providers,
        corrections.providers,
    )
)
del importlib
del itertools
