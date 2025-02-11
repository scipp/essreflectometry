# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc

from ..reflectometry.conversions import reflectometry_q
from ..reflectometry.types import (
    BeamDivergenceLimits,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)
from .types import CoordTransformationGraph


def theta(position, sample_position, sample_rotation, incident_angle_of_center_of_beam):
    '''
    Angle of reflection.

    Computes the angle between the scattering direction of
    the neutron and the sample surface.
    '''
    p = position - sample_position.to(unit=position.unit)
    return (
        sc.atan2(y=p.fields.x, x=p.fields.z)
        - sample_rotation.to(unit='rad')
        - incident_angle_of_center_of_beam.to(unit='rad')
    )


def divergence_angle(
    theta,
    sample_rotation,
):
    """
    Angle between the scattering direction and
    the ray from the sample to the center of the detector.
    """
    return theta - sample_rotation.to(unit=theta.unit)


def wavelength(
    event_time_offset,
    # Other inputs
):
    "Converts event_time_offset to wavelength"
    # Use frame unwrapping from scippneutron
    pass


def coordinate_transformation_graph() -> CoordTransformationGraph:
    return {
        "wavelength": wavelength,
        "theta": theta,
        "divergence_angle": divergence_angle,
        "Q": reflectometry_q,
        "L1": lambda source_position, sample_position: sc.norm(
            sample_position - source_position
        ),  # + extra correction for guides?
        "L2": lambda position, sample_position: sc.norm(position - sample_position),
        "incident_angle_of_center_of_beam": lambda: sc.scalar(1.7, unit='deg').to(
            unit='rad'
        ),
    }


def add_coords(
    da: sc.DataArray,
    graph: dict,
) -> sc.DataArray:
    "Adds scattering coordinates to the raw detector data."
    return da.transform_coords(
        ("wavelength", "theta", "divergence_angle", "Q", "L1", "L2"),
        graph,
        rename_dims=False,
        keep_intermediate=False,
        keep_aliases=False,
    )


def _not_between(v, a, b):
    return (v < a) | (v > b)


def add_masks(
    da: sc.DataArray,
    ylim: YIndexLimits,
    zlims: ZIndexLimits,
    bdlim: BeamDivergenceLimits,
    wbins: WavelengthBins,
) -> sc.DataArray:
    """
    Masks the data by ranges in the detector
    coordinates ``z`` and ``y``, and by the divergence of the beam,
    and by wavelength.
    """
    da = da.assign_masks(
        stripe_range=_not_between(da.coords["stripe"], *ylim),
        z_range=_not_between(da.coords["z_index"], *zlims),
        divergence_too_large=_not_between(
            da.coords["divergence_angle"],
            bdlim[0].to(unit=da.coords["divergence_angle"].unit),
            bdlim[1].to(unit=da.coords["divergence_angle"].unit),
        ),
    )
    da = da.bins.assign_masks(
        wavelength=_not_between(
            da.bins.coords['wavelength'],
            wbins[0],
            wbins[-1],
        ),
    )
    return da


providers = (coordinate_transformation_graph,)
