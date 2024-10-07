# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)

from .types import ReflectivityOverQ, ReflectivityOverQ1D


def average_over_detectors(iofq: ReflectivityOverQ) -> ReflectivityOverQ1D:
    """Average over all detector pixels."""
    return iofq.mean(dim='detector_number')


providers = (average_over_detectors,)
