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
    c = sc.constants.g * sc.constants.m_n**2 / sc.constants.h**2
    out = (c * L2 * wavelength**2).to(unit='dimensionless') + sc.sin(
        divergence_angle + detector_rotation
    )
    out = sc.asin(out, out=out)
    out -= sample_rotation
    return out


def angle_of_divergence(theta, sample_rotation, angle_to_center_of_beam):
    return (
        theta
        - sample_rotation
        - angle_to_center_of_beam
        - sc.scalar(0.245, unit='deg').to(unit='rad')
    )


def wavelength(
    event_time_offset, divergence_angle, L1, L2, chopper_phase, chopper_frequency
):
    out = event_time_offset.to(unit="ns", dtype="float64", copy=True)
    unit = out.bins.unit
    tau = (1 / (2 * chopper_frequency)).to(unit=unit)
    tof_offset = tau * chopper_phase / (180.0 * sc.units.deg)

    minimum = -tof_offset
    frame_bound = tau - tof_offset
    maximum = 2 * tau - tof_offset

    out += sc.where(
        (minimum < event_time_offset) & (event_time_offset < frame_bound),
        tof_offset,
        sc.where(
            (frame_bound < event_time_offset) & (event_time_offset < maximum),
            tof_offset - tau,
            sc.scalar(np.nan, unit=unit),
        ),
    )
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
