# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import numpy as np
import scipp as sc

from ..reflectometry.conversions import reflectometry_q
from ..reflectometry.corrections import footprint_on_sample
from ..reflectometry.types import (
    BeamDivergenceLimits,
    BeamSize,
    RawDetectorData,
    ReducibleData,
    RunType,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)


def theta(wavelength, divergence_angle, L2, sample_rotation, detector_rotation):
    '''
    Angle of reflection.

    Computes the angle between the scattering direction of
    the neutron and the sample surface.

    :math:`\\gamma^*` denotes the angle between the scattering direction
    and the horizontal plane.
    :math:`\\gamma` denotes the angle between the ray from sample position
    to detection position
    and the horizontal plane.
    :math:`L_2` is the length of the ray from sample position to detector position.
    :math:`v` is the velocity of the neutron at the sample.
    :math:`t` is the travel time from sample to detector.

    The parabolic trajectory of the neutron satisfies

    .. math::

        \\sin(\\gamma) L_2 = \\sin(\\gamma^*) v t - \\frac{g}{2} t^2

    and

    .. math::

        \\cos(\\gamma) L_2 = \\cos(\\gamma^*) vt

    where :math:`g` is the gravitational acceleration.

    The second equation tells us that the approximation :math:`L_2=vt`
    will have a small error if :math:`\\gamma` is close to 0 and
    the difference between :math:`\\gamma` and :math:`\\gamma^*` is small.

    Using this approximation we can solve the first equation,
    and by expressing :math:`v` in terms of the wavelength we get

    .. math::

        \\sin(\\gamma^*) =
        \\sin(\\gamma) + \\frac{g}{2} \\frac{L_2 \\lambda^2 h^2}{m_n^2}.

    Finally, the scattering angle is obtained by subtracting the sample rotation
    relative to the horizontal plane.
    '''
    c = sc.constants.g * sc.constants.m_n**2 / sc.constants.h**2
    out = (c * L2 * wavelength**2).to(unit='dimensionless') + sc.sin(
        divergence_angle + detector_rotation
    )
    out = sc.asin(out, out=out)
    out -= sample_rotation
    return out


def angle_of_divergence(theta, sample_rotation, angle_to_center_of_beam):
    """
    Difference between the incident angle and the center of the incident beam.
    Useful for filtering parts of the beam that have too high divergence.

    On the Amor instrument this is always in the interval [-0.75 deg, 0.75 deg],
    but the divergence of the incident beam can be made lower.
    """
    return theta - sample_rotation - angle_to_center_of_beam


def wavelength(
    event_time_offset, divergence_angle, L1, L2, chopper_phase, chopper_frequency
):
    "Converts event_time_offset to wavelength using the chopper settings."
    out = event_time_offset.to(unit="ns", dtype="float64", copy=True)
    unit = out.bins.unit
    tau = (1 / (2 * chopper_frequency)).to(unit=unit)
    tof_offset = tau * chopper_phase / (180.0 * sc.units.deg)

    minimum = -tof_offset
    frame_bound = tau - tof_offset
    maximum = 2 * tau - tof_offset

    # Frame unwrapping
    out += sc.where(
        (minimum < event_time_offset) & (event_time_offset < frame_bound),
        tof_offset,
        sc.where(
            (frame_bound < event_time_offset) & (event_time_offset < maximum),
            tof_offset - tau,
            sc.scalar(np.nan, unit=unit),
        ),
    )
    # Correction for path length through guides being different
    # depending on incident angle.
    out -= (divergence_angle.to(unit="deg") / (180.0 * sc.units.deg)) * tau
    out *= (sc.constants.h / sc.constants.m_n) / (L1 + L2)
    return out.to(unit='angstrom', copy=False)


def _not_between(v, a, b):
    return (v < a) | (v > b)


def add_common_coords_and_masks(
    da: RawDetectorData[RunType],
    ylim: YIndexLimits,
    zlims: ZIndexLimits,
    bdlim: BeamDivergenceLimits,
    wbins: WavelengthBins,
    beam_size: BeamSize[RunType],
) -> ReducibleData[RunType]:
    "Adds coords and masks that are useful for both reference and sample measurements."
    da = da.transform_coords(
        ("wavelength", "theta", "angle_of_divergence", "Q"),
        {
            "divergence_angle": "pixel_divergence_angle",
            "wavelength": wavelength,
            "theta": theta,
            "angle_of_divergence": angle_of_divergence,
            "Q": reflectometry_q,
        },
        rename_dims=False,
        keep_intermediate=False,
    )
    da.masks["stripe_range"] = _not_between(da.coords["stripe"], *ylim)
    da.masks['z_range'] = _not_between(da.coords["z_index"], *zlims)
    da.bins.masks["divergence_too_large"] = _not_between(
        da.bins.coords["angle_of_divergence"],
        bdlim[0].to(unit=da.bins.coords["angle_of_divergence"].bins.unit),
        bdlim[1].to(unit=da.bins.coords["angle_of_divergence"].bins.unit),
    )
    da.bins.masks['wavelength'] = _not_between(
        da.bins.coords['wavelength'],
        wbins[0],
        wbins[-1],
    )
    da /= footprint_on_sample(
        da.bins.coords['theta'],
        beam_size,
        da.coords['sample_size'],
    )
    return da


providers = (add_common_coords_and_masks,)
