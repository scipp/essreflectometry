# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
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


def stitch_reflecivity_curves(curves, qgrid):
    '''Stitches the curves by scaling each except the first by a factor.
    The scaling factors are determined by a maximum likelihood estimate
    (assuming the errors are normal distributed).
    '''

    def evaluate_on_qgrid(f):
        return np.stack(
            [
                np.interp(
                    qgrid.values,
                    sc.midpoints(c.coords['Q']).values,
                    f(c),
                    left=np.nan,
                    right=np.nan,
                )
                for c in curves
            ]
        )

    r = evaluate_on_qgrid(lambda c: c.data.values)
    v = evaluate_on_qgrid(lambda c: c.data.variances)

    def cost(ps):
        ps = np.concatenate([[1.0], ps])
        rs = ps[:, None] * r
        ss = ps[:, None] ** 2 * v
        ss[ss == 0] = np.nan
        iss = 1 / ss
        m = np.nansum(rs * iss, axis=0) / np.nansum(iss, axis=0)
        return np.nansum((rs - m) ** 2 * iss)

    sol = opt.minimize(cost, [1.0] * (len(curves) - 1))
    return [p * c for p, c in zip((1.0, *sol.x), curves, strict=True)]


def combine_curves(curves, qgrid, how='mean'):
    '''Combines the given curves by interpolating them
    on a grid and merging them by the requested method.
    The default method is a weighted mean where the weights
    are proportional to the variances.'''

    def evaluate_on_qgrid(f):
        return np.stack(
            [
                np.interp(
                    sc.midpoints(qgrid).values,
                    sc.midpoints(c.coords['Q']).values,
                    f(c),
                    left=np.nan,
                    right=np.nan,
                )
                for c in curves
            ]
        )

    r = evaluate_on_qgrid(lambda c: c.data.values)
    v = evaluate_on_qgrid(lambda c: c.data.variances)

    if how == 'mean':
        v[v == 0] = np.nan
        iv = 1.0 / v
        m = np.nansum(r * iv, axis=0) / np.nansum(iv, axis=0)
        mv = 1 / np.nansum(iv, axis=0)
        return sc.DataArray(
            data=sc.array(dims='Q', values=m, variances=mv), coords={'Q': qgrid}
        )
    return NotImplementedError
