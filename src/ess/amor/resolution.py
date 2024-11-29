# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc

from ..reflectometry.types import (
    DetectorSpatialResolution,
    QBins,
    QResolution,
    Reference,
    SampleRun,
    SampleSize,
)
from .types import (
    AngularResolution,
    ChopperSeparation,
    SampleSizeResolution,
    WavelengthResolution,
)

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


def wavelength_resolution(
    L1,
    L2,
    chopper_separation,
):
    """
    Find the wavelength resolution contribution as described in Section 4.3.3 of the
    Amor publication (doi: 10.1016/j.nima.2016.03.007).

    Parameters
    ----------

    L1:
        Distance from midpoint between choppers to sample.
    L2:
        Distance from sample to detector.
    chopper_separation:
        Distance between choppers.

    Returns
    -------
    :
        The angular resolution variable, as standard deviation.
    """
    return fwhm_to_std(chopper_separation / (L1 + L2))


def _wavelength_resolution(
    da: Reference, chopper_separation: ChopperSeparation[SampleRun]
) -> WavelengthResolution:
    return wavelength_resolution(
        L1=da.coords['L1'], L2=da.coords['L2'], chopper_separation=chopper_separation
    )


def sample_size_resolution(
    L2,
    sample_size,
):
    """
    The resolution from the projected sample size, where it may be bigger
    than the detector pixel resolution as described in Section 4.3.3 of the Amor
    publication (doi: 10.1016/j.nima.2016.03.007).

    Parameters
    ----------
    L2:
        Distance from sample to detector.
    sample_size:
        Size of sample.

    Returns
    -------
    :
        Standard deviation of contribution from the sample size.
    """
    return fwhm_to_std(sample_size / L2.to(unit=sample_size.unit))


def _sample_size_resolution(
    da: Reference, sample_size: SampleSize[SampleRun]
) -> SampleSizeResolution:
    return sample_size_resolution(L2=da.coords['L2'], sample_size=sample_size)


def angular_resolution(
    theta,
    L2,
    detector_spatial_resolution,
):
    """
    Determine the angular resolution as described in Section 4.3.3 of the Amor
    publication (doi: 10.1016/j.nima.2016.03.007).

    Parameters
    ----------
    theta:
        Angle of reflection.
    L2:
        Distance between sample and detector.
    detector_spatial_resolution:
        FWHM of detector pixel resolution.

    Returns
    -------
    :
        Angular resolution standard deviation
    """
    return (
        fwhm_to_std(
            sc.atan(
                detector_spatial_resolution
                / L2.to(unit=detector_spatial_resolution.unit)
            )
        ).to(unit=theta.unit)
        / theta
    )


def _angular_resolution(
    da: Reference, detector_spatial_resolution: DetectorSpatialResolution[SampleRun]
) -> AngularResolution:
    return angular_resolution(
        theta=da.coords['theta'],
        L2=da.coords['L2'],
        detector_spatial_resolution=detector_spatial_resolution,
    )


def sigma_Q(
    ref: Reference,
    angular_resolution: AngularResolution,
    wavelength_resolution: WavelengthResolution,
    sample_size_resolution: SampleSizeResolution,
    qbins: QBins,
) -> QResolution:
    """
    Combine all of the components of the resolution and add Q contribution.

    Parameters
    ----------
    angular_resolution:
        Angular resolution contribution.
    wavelength_resolution:
        Wavelength resolution contribution.
    sample_size_resolution:
        Sample size resolution contribution.
    q_bins:
        Q-bin values.

    Returns
    -------
    :
        Combined resolution function.
    """
    h = ref.flatten(to='Q').hist(Q=qbins)
    s = (
        (
            ref
            * (
                angular_resolution**2
                + wavelength_resolution**2
                + sample_size_resolution**2
            )
            * ref.coords['Q'] ** 2
        )
        .flatten(to='Q')
        .hist(Q=qbins)
    )
    return sc.sqrt(sc.values(s / h))


providers = (
    sigma_Q,
    _angular_resolution,
    _wavelength_resolution,
    _sample_size_resolution,
)
