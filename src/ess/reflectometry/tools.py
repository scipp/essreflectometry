# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
from collections.abc import Sequence
from itertools import chain

import numpy as np
import scipp as sc
import scipy.optimize as opt

_STD_TO_FWHM = sc.scalar(2.0) * sc.sqrt(sc.scalar(2.0) * sc.log(sc.scalar(2.0)))


def fwhm_to_std(fwhm: sc.Variable) -> sc.Variable:
    """
    Convert from full-width half maximum to standard deviation.

    Parameters
    ----------
    fwhm:
        Full-width half maximum.

    Returns
    -------
    :
        Standard deviation.
    """
    # Enables the conversion from full width half
    # maximum to standard deviation
    return fwhm / _STD_TO_FWHM


def std_to_fwhm(std: sc.Variable) -> sc.Variable:
    """
    Convert from standard deviation to full-width half maximum.

    Parameters
    ----------
    std:
        Standard deviation.

    Returns
    -------
    :
        Full-width half maximum.
    """
    # Enables the conversion from full width half
    # maximum to standard deviation
    return std * _STD_TO_FWHM


def linlogspace(
    dim: str,
    edges: list | np.ndarray,
    scale: list | str,
    num: list | int,
    unit: str | None = None,
) -> sc.Variable:
    """
    Generate a 1d array of bin edges with a mixture of linear and/or logarithmic
    spacings.

    Examples:

    - Create linearly spaced edges (equivalent to `sc.linspace`):
        linlogspace(dim='x', edges=[0.008, 0.08], scale='linear', num=50, unit='m')
    - Create logarithmically spaced edges (equivalent to `sc.geomspace`):
        linlogspace(dim='x', edges=[0.008, 0.08], scale='log', num=50, unit='m')
    - Create edges with a linear and a logarithmic part:
        linlogspace(dim='x', edges=[1, 3, 8], scale=['linear', 'log'], num=[16, 20])

    Parameters
    ----------
    dim:
        The dimension of the output Variable.
    edges:
        The edges for the different parts of the mesh.
    scale:
        A string or list of strings specifying the scaling for the different
        parts of the mesh. Possible values for the scaling are `"linear"` and `"log"`.
        If a list is supplied, the length of the list must be one less than the length
        of the `edges` parameter.
    num:
        An integer or a list of integers specifying the number of points to use
        in each part of the mesh. If a list is supplied, the length of the list must be
        one less than the length of the `edges` parameter.
    unit:
        The unit of the output Variable.

    Returns
    -------
    :
        Lin-log spaced Q-bin edges.
    """
    if not isinstance(scale, list):
        scale = [scale]
    if not isinstance(num, list):
        num = [num]
    if len(scale) != len(edges) - 1:
        raise ValueError(
            "Sizes do not match. The length of edges should be one "
            "greater than scale."
        )

    funcs = {"linear": sc.linspace, "log": sc.geomspace}
    grids = []
    for i in range(len(edges) - 1):
        # Skip the leading edge in the piece when concatenating
        start = int(i > 0)
        mesh = funcs[scale[i]](
            dim=dim, start=edges[i], stop=edges[i + 1], num=num[i] + start, unit=unit
        )
        grids.append(mesh[dim, start:])

    return sc.concat(grids, dim)


def _sort_by(a, by):
    return [x for x, _ in sorted(zip(a, by, strict=True), key=lambda x: x[1])]


def _find_interval_overlaps(intervals):
    '''Returns the intervals where at least
    two or more of the provided intervals
    are overlapping.'''
    edges = list(chain.from_iterable(intervals))
    is_start_edge = list(chain.from_iterable((True, False) for _ in intervals))
    edges_sorted = sorted(edges)
    is_start_edge_sorted = _sort_by(is_start_edge, edges)

    number_overlapping = 0
    overlap_intervals = []
    for x, is_start in zip(edges_sorted, is_start_edge_sorted, strict=True):
        if number_overlapping == 1 and is_start:
            start = x
        if number_overlapping == 2 and not is_start:
            overlap_intervals.append((start, x))
        if is_start:
            number_overlapping += 1
        else:
            number_overlapping -= 1
    return overlap_intervals


def _searchsorted(a, v):
    for i, e in enumerate(a):
        if e > v:
            return i
    return len(a)


def _create_qgrid_where_overlapping(qgrids):
    '''Given a number of Q-grids, construct a new grid
    covering the regions where (any two of the) provided grids overlap.'''
    pieces = []
    for start, end in _find_interval_overlaps([(q.min(), q.max()) for q in qgrids]):
        interval_sliced_from_qgrids = [
            q[max(_searchsorted(q, start) - 1, 0) : _searchsorted(q, end) + 1]
            for q in qgrids
        ]
        densest_grid_in_interval = max(interval_sliced_from_qgrids, key=len)
        pieces.append(densest_grid_in_interval)
    return sc.concat(pieces, dim='Q')


def _interpolate_on_qgrid(curves, grid):
    return sc.concat(
        [sc.lookup(c, grid.dim)[sc.midpoints(grid)] for c in curves], dim='curves'
    )


def scale_reflectivity_curves_to_overlap(
    curves: Sequence[sc.DataArray],
    return_scaling_factors=False,
) -> list[sc.DataArray] | list[sc.scalar]:
    '''Make the curves overlap by scaling all except the first by a factor.
    The scaling factors are determined by a maximum likelihood estimate
    (assuming the errors are normal distributed).

    All curves must be have the same unit for data and the Q-coordinate.

    Parameters
    ---------
    curves:
        the reflectivity curves that should be scaled together
    return_scaling_factor:
        If True the return value of the function
        is a list of the scaling factors that should be applied.
        If False (default) the function returns the scaled curves.

    Returns
    ---------
    :
        A list of scaled reflectivity curves or a list of scaling factors.
    '''
    if len({c.data.unit for c in curves}) != 1:
        raise ValueError('The reflectivity curves must have the same unit')
    if len({c.coords['Q'].unit for c in curves}) != 1:
        raise ValueError('The Q-coordinates must have the same unit for each curve')

    qgrid = _create_qgrid_where_overlapping([c.coords['Q'] for c in curves])

    r = _interpolate_on_qgrid(map(sc.values, curves), qgrid).values
    v = _interpolate_on_qgrid(map(sc.variances, curves), qgrid).values

    def cost(scaling_factors):
        scaling_factors = np.concatenate([[1.0], scaling_factors])[:, None]
        r_scaled = scaling_factors * r
        v_scaled = scaling_factors**2 * v
        v_scaled[v_scaled == 0] = np.nan
        inv_v_scaled = 1 / v_scaled
        r_avg = np.nansum(r_scaled * inv_v_scaled, axis=0) / np.nansum(
            inv_v_scaled, axis=0
        )
        return np.nansum((r_scaled - r_avg) ** 2 * inv_v_scaled)

    sol = opt.minimize(cost, [1.0] * (len(curves) - 1))
    scaling_factors = (1.0, *sol.x)
    if return_scaling_factors:
        return [sc.scalar(x) for x in scaling_factors]
    return [
        scaling_factor * curve
        for scaling_factor, curve in zip(scaling_factors, curves, strict=True)
    ]


def combine_curves(
    curves: Sequence[sc.DataArray],
    qgrid: sc.Variable | None = None,
) -> sc.DataArray:
    '''Combines the given curves by interpolating them
    on a grid and merging them by the requested method.
    The default method is a weighted mean where the weights
    are proportional to the variances.

    Unless the curves are already scaled correctly they might
    need to be scaled using :func:`scale_reflectivity_curves_to_overlap`.

    All curves must be have the same unit for data and the Q-coordinate.

    Parameters
    ----------
    curves:
        the reflectivity curves that should be combined
    qgrid:
        the Q-grid of the resulting combined reflectivity curve

    Returns
    ---------
    :
        A data array representing the combined reflectivity curve
    '''
    if len({c.data.unit for c in curves}) != 1:
        raise ValueError('The reflectivity curves must have the same unit')
    if len({c.coords['Q'].unit for c in curves}) != 1:
        raise ValueError('The Q-coordinates must have the same unit for each curve')

    r = _interpolate_on_qgrid(map(sc.values, curves), qgrid).values
    v = _interpolate_on_qgrid(map(sc.variances, curves), qgrid).values

    v[v == 0] = np.nan
    inv_v = 1.0 / v
    r_avg = np.nansum(r * inv_v, axis=0) / np.nansum(inv_v, axis=0)
    v_avg = 1 / np.nansum(inv_v, axis=0)
    return sc.DataArray(
        data=sc.array(
            dims='Q',
            values=r_avg,
            variances=v_avg,
            unit=next(iter(curves)).data.unit,
        ),
        coords={'Q': qgrid},
    )
