# Copyright (c) 2025 Scipp contributors (https://github.com/scipp)
from ..reflectometry.conversions import reflectometry_q
from ..reflectometry.types import (
    QResolution,
    ReducibleData,
    Reference,
    ReferenceRun,
    Sample,
    SampleRun,
)


def evaluate_reference(
    reference: ReducibleData[ReferenceRun],
    sample: ReducibleData[SampleRun],
    qresolution: QResolution,
) -> Reference:
    """
    Adds a :math:`Q`. The coordinate is computed as if the data came from
    the sample measurement, that is, they use the ``sample_rotation``
    of the sample measurement.
    """
    ref = reference.copy()
    ref.coords.pop("theta")
    ref.bins.coords['Q'] = reflectometry_q(
        wavelength=ref.bins.coords['wavelength'], theta=sample.coords['theta']
    )
    ref.bins.coords['Q_resolution'] = qresolution * ref.bins.coords['Q']
    return ref.bins.concat()


def evaluate_sample(
    reference: ReducibleData[ReferenceRun],
    sample: ReducibleData[SampleRun],
) -> Sample:
    """
    Adds a :math:`Q`.
    """
    sample = sample.copy()
    sample.bins.coords['Q'] = reflectometry_q(
        wavelength=sample.bins.coords['wavelength'], theta=sample.coords['theta']
    )
    return sample.bins.concat()


providers = (
    evaluate_reference,
    evaluate_sample,
)
