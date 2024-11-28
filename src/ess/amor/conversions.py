# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import numpy as np
import scipp as sc

from ..reflectometry.types import (
    BeamDivergenceLimits,
    RawDetectorData,
    ReducibleData,
    RunType,
    YIndexLimits,
)


def angle_of_reflection(
    wavelength, divergence_angle, L2, sample_rotation, detector_rotation
):
    c = sc.constants.g * sc.constants.m_n**2 / sc.constants.h**2
    out = (c * L2 * wavelength**2).to(unit='dimensionless') + sc.sin(
        divergence_angle + detector_rotation
    )
    out = sc.asin(out, out=out)
    out -= sample_rotation
    return out


def angle_of_divergence(angle_of_reflection, sample_rotation):
    return angle_of_reflection - sample_rotation


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


def add_common_coords_and_masks(
    da: RawDetectorData[RunType],
    ylim: YIndexLimits,
    bdlim: BeamDivergenceLimits,
) -> ReducibleData[RunType]:
    da.bins.coords["wavelength"] = wavelength(
        event_time_offset=da.bins.coords["event_time_offset"],
        divergence_angle=da.coords["pixel_divergence_angle"],
        L1=da.coords["L1"],
        L2=da.coords["L2"],
        chopper_phase=da.coords["chopper_phase"],
        chopper_frequency=da.coords["chopper_frequency"],
    )
    da.bins.coords["angle_of_reflection"] = angle_of_reflection(
        wavelength=da.bins.coords["wavelength"],
        divergence_angle=da.coords["pixel_divergence_angle"],
        L2=da.coords["L2"],
        sample_rotation=da.coords["sample_rotation"],
        detector_rotation=da.coords["detector_rotation"],
    )
    da.bins.coords["angle_of_divergence"] = angle_of_divergence(
        angle_of_reflection=da.bins.coords["angle_of_reflection"],
        sample_rotation=da.coords["sample_rotation"],
    )
    da.masks["stripe_range"] = (da.coords["stripe"] < ylim[0]) | (
        da.coords["stripe"] > ylim[1]
    )
    da.bins.masks["divergence_too_large"] = (
        da.bins.coords["angle_of_divergence"]
        < bdlim[0].to(unit=da.bins.coords["angle_of_divergence"].bins.unit)
    ) | (
        da.bins.coords["angle_of_divergence"]
        > bdlim[1].to(unit=da.bins.coords["angle_of_divergence"].bins.unit)
    )
    # Footprint correction
    da /= sc.sin(da.bins.coords['angle_of_reflection'])
    return da


providers = (add_common_coords_and_masks,)
