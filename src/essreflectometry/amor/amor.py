# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)

import scipp as sc

from .. import reflectometry
from . import calibrations, conversions, data
from . import normalize as normalize_amor
from . import resolution
from .beamline import make_beamline
from .load import load as load_amor
from .types import (
    CalibratedReference,
    CountsByMomentumTransfer,
    Experiment,
    FootprintCorrected,
    MomentumTransferBins,
    Normalized,
    Raw,
    Reference,
    Rotation,
    Run,
    Sample,
    WavelengthBins,
    WithQResolution,
)


def load(filename: Run, rotation: Rotation[Run]) -> Raw[Run]:
    '''Load raw experiment data + geometry'''
    da = load_amor(
        data.get_path(filename),
        orso=None,
        beamline=make_beamline(sample_rotation=rotation),
    )
    da.coords['position'].fields.y += da.coords['position'].fields.z * sc.tan(
        2.0 * da.coords['sample_rotation'] - (0.955 * sc.units.deg)
    )
    return da


def compute_coordinates_based_on_time_of_flight(
    da: Raw[Run],
    wavelength_bins: WavelengthBins,
    momentum_transfer_bins: MomentumTransferBins,
) -> Experiment[Run]:
    graph = conversions.specular_reflection()
    da = reflectometry.conversions.tof_to_wavelength(
        da,
        wavelength_bins,
        graph=graph,
    )
    da = reflectometry.conversions.wavelength_to_theta(
        da,
        graph=graph,
    )
    return reflectometry.conversions.theta_to_q(
        da,
        q_edges=momentum_transfer_bins,
        graph=graph,
    )


def correct_for_footprint(da: Experiment[Run]) -> FootprintCorrected[Run]:
    '''Correct counts by footprint on detector'''
    return reflectometry.corrections.footprint_correction(da)


def calibrate_reference(da: FootprintCorrected[Reference]) -> CalibratedReference:
    '''Correct counts by footprint on detector'''
    reference_q_summed = reflectometry.conversions.sum_bins(da)
    return calibrations.supermirror_calibration(reference_q_summed)


def add_momentum_transfer_resolution(da: FootprintCorrected[Sample]) -> WithQResolution:
    da.coords['wavelength_resolution'] = resolution.wavelength_resolution(
        chopper_1_position=da.coords['source_chopper_1'].value['position'],
        chopper_2_position=da.coords['source_chopper_2'].value['position'],
        pixel_position=da.coords['position'],
    )
    da.bins.coords['angular_resolution'] = resolution.angular_resolution(
        pixel_position=da.coords['position'],
        theta=da.bins.coords['theta'],
        detector_spatial_resolution=da.coords['detector_spatial_resolution'],
    )
    da.coords['sample_size_resolution'] = resolution.sample_size_resolution(
        pixel_position=da.coords['position'], sample_size=da.coords['sample_size']
    )
    return da


def normalize(da: WithQResolution, ref: CalibratedReference) -> Normalized:
    '''Normalize counts by calibrated reference'''
    summed = reflectometry.conversions.sum_bins(da)
    norm = reflectometry.corrections.normalize_by_counts(summed)
    reference_norm = reflectometry.corrections.normalize_by_counts(ref)
    da = normalize_amor.normalize_by_supermirror(norm, reference_norm)
    da.coords['sigma_Q'] = resolution.sigma_Q(
        angular_resolution=da.coords['angular_resolution'],
        wavelength_resolution=da.coords['wavelength_resolution'],
        sample_size_resolution=da.coords['sample_size_resolution'],
        q_bins=da.coords['Q'],
    )
    return da


def histogram_over_momentum_transfer(da: Normalized) -> CountsByMomentumTransfer:
    return da.mean('detector_number')


providers = [
    load,
    compute_coordinates_based_on_time_of_flight,
    correct_for_footprint,
    calibrate_reference,
    normalize,
    add_momentum_transfer_resolution,
]
