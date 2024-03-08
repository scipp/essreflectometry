# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Scipp contributors (https://github.com/scipp)
"""ORSO utilities for Amor."""
from typing import Optional

import numpy as np
import scipp as sc
from orsopy.fileio import base as orso_base
from orsopy.fileio import data_source as orso_data_source
from orsopy.fileio.orso import Column, Orso, OrsoDataset

from ..reflectometry.orso import (
    OrsoDataSource,
    OrsoInstrument,
    OrsoIofQDataset,
    OrsoReduction,
)
from ..reflectometry.types import (
    NormalizedIofQ1D,
    QResolution,
    Sample,
    ThetaData,
    WavelengthData,
)


def build_orso_instrument(
    events_in_wavelength: WavelengthData[Sample], events_in_theta: ThetaData[Sample]
) -> OrsoInstrument:
    """Build ORSO instrument metadata from intermediate reduction results for Amor.

    This assumes specular reflection and sets the incident angle equal to the computed
    scattering angle.
    """
    return OrsoInstrument(
        orso_data_source.InstrumentSettings(
            wavelength=orso_base.ValueRange(
                *_limits_of_coord(events_in_wavelength, 'wavelength')
            ),
            incident_angle=orso_base.ValueRange(
                *_limits_of_coord(events_in_theta, 'theta')
            ),
            polarization=None,  # TODO how can we determine this from the inputs?
        )
    )


def build_orso_iofq_dataset(
    iofq: NormalizedIofQ1D,
    sigma_q: QResolution,
    data_source: OrsoDataSource,
    reduction: OrsoReduction,
) -> OrsoIofQDataset:
    """Build an ORSO dataset for reduced I-of-Q data and associated metadata."""
    header = Orso(
        data_source=data_source,
        reduction=reduction,
        columns=[
            Column('Qz', '1/angstrom', 'wavevector transfer'),
            Column('R', None, 'reflectivity'),
            Column('sR', None, 'standard deviation of reflectivity'),
            Column(
                'sQz',
                '1/angstrom',
                'standard deviation of wavevector transfer resolution',
            ),
        ],
    )

    qz = iofq.coords['Q'].to(unit='1/angstrom', copy=False)
    if iofq.coords.is_edges('Q'):
        qz = sc.midpoints(qz)
    r = sc.values(iofq.data)
    sr = sc.stddevs(iofq.data)
    sqz = sigma_q.to(unit='1/angstrom', copy=False)
    data = (qz, r, sr, sqz)

    return OrsoIofQDataset(
        OrsoDataset(header, np.column_stack([_extract_values_array(d) for d in data]))
    )


def _extract_values_array(var: sc.Variable) -> np.ndarray:
    if var.variances is not None:
        raise sc.VariancesError(
            "ORT columns must not have variances. "
            "Store the uncertainties as standard deviations in a separate column."
        )
    if var.ndim != 1:
        raise sc.DimensionError(f"ORT columns must be one-dimensional, got {var.sizes}")
    return var.values


def _limits_of_coord(
    data: sc.DataArray, name: str
) -> Optional[tuple[float, float, str]]:
    if (coord := _get_coord(data, name)) is None:
        return None
    min_ = coord.min().value
    max_ = coord.max().value
    # Explicit conversions to float because orsopy does not like np.float* types.
    return float(min_), float(max_), _ascii_unit(min_)


def _get_coord(data: sc.DataArray, name: str) -> Optional[sc.Variable]:
    try:
        return data.coords[name]
    except KeyError:
        try:
            return data.bins.coords[name]
        except (KeyError, TypeError):
            # KeyError if coord does not exist, TypeError if data is not binned.
            return None


def _ascii_unit(unit: sc.Unit) -> str:
    unit = str(unit)
    if unit == 'Å':
        return 'angstrom'
    return unit


providers = (build_orso_instrument, build_orso_iofq_dataset)
