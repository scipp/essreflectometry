# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc

# from ..reflectometry import orso
from ..types import CalibratedReference, Histogrammed, Reference
from .types import Alpha, CriticalEdge, MValue, Supermirror


def supermirror_calibration(
    data_array: Histogrammed[Reference],
    m_value: Supermirror[MValue],
    critical_edge: Supermirror[CriticalEdge],
    alpha: Supermirror[Alpha],
) -> CalibratedReference:
    """
    Calibrate supermirror measurements

    Parameters
    ----------
    data_array:
        Data array to get q-bins/values from.
    m_value:
        m-value for the supermirror.
    critical_edge:
        Supermirror critical edge.
    alpha:
        Supermirror alpha value.
    Returns
    -------
    :
        Calibrated supermirror measurement.
    """
    calibration = calibration_factor(data_array, m_value, critical_edge, alpha)
    data_array_cal = data_array * calibration
    # TODO
    # try:
    #    data_array_cal.attrs['orso'].value.reduction.corrections += [
    #        'supermirror calibration'
    #    ]
    # except KeyError:
    #    orso.not_found_warning()
    return Histogrammed[CalibratedReference](data_array_cal)


def calibration_factor(
    data_array: sc.DataArray,
    m_value: sc.Variable,
    critical_edge: sc.Variable,
    alpha: sc.Variable,
) -> sc.Variable:
    """
    Return the calibration factor for the supermirror.

    Parameters
    ----------
    data_array:
        Data array to get q-bins/values from.
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
    q = data_array.coords['Q']
    if data_array.coords.is_edges('Q'):
        q = sc.midpoints(q)
    max_q = m_value * critical_edge
    lim = (q < critical_edge).astype(float)
    lim.unit = 'one'
    nq = 1.0 / (1.0 - alpha * (q - critical_edge))
    calibration_factor = sc.where(q < max_q, lim + (1 - lim) * nq, sc.scalar(1.0))
    return calibration_factor


providers = [supermirror_calibration]
