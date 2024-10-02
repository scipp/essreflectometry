# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)

import scipp as sc

from .supermirror import SupermirrorReflectivityCorrection
from .tools import fwhm_to_std
from .types import (
    BeamSize,
    CorrectedData,
    EventData,
    FootprintCorrection,
    IdealReferenceIntensity,
    MaskedData,
    ProtonCurrent,
    ProtonCurrentCorrection,
    RawEventData,
    ReferenceIntensity,
    ReferenceRun,
    RunType,
    SampleSize,
    WavelengthBins,
)


def footprint_correction(
    data_array: MaskedData[RunType],
    beam_size: BeamSize[RunType],
    sample_size: SampleSize[RunType],
) -> FootprintCorrection[RunType]:
    """
    Corrects the event weights by the fraction of the beam hitting the sample.
    Depends on :math:`\\theta`.

    Parameters
    ----------
    data_array:
        Data array to perform footprint correction on.
    beam_size:
        Full width half maximum of the beam.
    sample_size:
        Size of the sample.
        TODO: check what sample size actually means. Is it the sample diameter? etc.

    Returns
    -------
    :
       Footprint corrected data array.
    """
    size_of_beam_on_sample = beam_size / sc.sin(data_array.bins.coords["theta"])
    footprint_scale = sc.erf(fwhm_to_std(sample_size / size_of_beam_on_sample))
    return FootprintCorrection[RunType](1 / footprint_scale)


def compute_reference_intensity(
    da: CorrectedData[ReferenceRun], wb: WavelengthBins
) -> ReferenceIntensity:
    """Creates a reference intensity map over (z_index, wavelength).
    Rationale:
        The intensity expressed in those variables should not vary
        with the experiment parameters (such as sample rotation).
        Therefore it can be used to normalize sample measurements.
    """
    b = da.bin(wavelength=wb, dim=set(da.dims) - set(da.coords["z_index"].dims))
    h = b.hist()
    h.masks["too_few_events"] = sc.stddevs(h.data) >= sc.values(h.data) / 2
    # Add a Q coordinate to each bin, the Q is not completely unique in every bin,
    # but it is close enough.
    h.coords["Q"] = b.bins.coords["Q"].bins.mean()
    return ReferenceIntensity(h)


def calibrate_reference(
    da: ReferenceIntensity, cal: SupermirrorReflectivityCorrection
) -> IdealReferenceIntensity:
    """Calibrates the reference intensity by the
    inverse of the supermirror reflectivity"""
    return IdealReferenceIntensity(da * cal)


def proton_current_correction(
    pc: ProtonCurrent[RunType], da: RawEventData[RunType]
) -> ProtonCurrentCorrection[RunType]:
    pcl = sc.lookup(pc, dim='time', mode='previous', fill_value=sc.scalar(float('nan')))
    proton_current_at_pulses = pcl(da.coords['event_time_zero'])
    return ProtonCurrentCorrection[RunType](1 / proton_current_at_pulses)


def correct_detector_binned_events(
    da: MaskedData[RunType], fp: FootprintCorrection[RunType]
) -> CorrectedData[RunType]:
    return CorrectedData[RunType](da * fp)


def correct_pulse_binned_events(
    da: RawEventData[RunType], pcc: ProtonCurrentCorrection[RunType]
) -> EventData[RunType]:
    return EventData[RunType](da)
    # If the proton current varies by a factor more than cutoff from the median
    # it is considered a bad pulse.
    cutoff = 5
    da.masks['too_low_proton_current'] = pcc > cutoff * sc.nanmedian(pcc)
    da.masks['undefined_proton_current'] = sc.isnan(pcc)
    return EventData[RunType](da * pcc)


providers = (
    footprint_correction,
    calibrate_reference,
    compute_reference_intensity,
    proton_current_correction,
    correct_detector_binned_events,
    correct_pulse_binned_events,
)
