# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc

from .supermirror.types import Alpha, CriticalEdge, MValue, SupermirrorCalibrationFactor

# from ..reflectometry import orso
from .types import QBins


def calibration_factor(
    qbins: QBins,
    m_value: MValue,
    critical_edge: CriticalEdge,
    alpha: Alpha,
) -> SupermirrorCalibrationFactor:
    """
    Return the calibration factor for the supermirror.

    Parameters
    ----------
    qbins:
        edges of binning of Q.
    m_value:
        m-value for the supermirror.
    critical_edge:
        Supermirror critical edge.
    alpha:
        Supermirror alpha value.

    Returns
    -------
    :
        Calibration factor at the midpoint of each Q-bin.
    """
    q = sc.midpoints(qbins)
    max_q = m_value * critical_edge
    lim = (q < critical_edge).astype(float)
    lim.unit = 'one'
    nq = 1.0 / (1.0 - alpha * (q - critical_edge))
    calibration_factor = sc.where(q < max_q, lim + (1 - lim) * nq, sc.scalar(1.0))
    return calibration_factor


providers = [calibration_factor]
