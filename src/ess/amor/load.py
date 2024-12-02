# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc

from ess.reduce import nexus

from ..reflectometry.load import load_nx
from ..reflectometry.types import (
    DetectorRotation,
    Filename,
    LoadedNeXusDetector,
    NeXusDetectorName,
    RawDetectorData,
    RunType,
    SampleRotation,
    SampleSize,
)
from .geometry import Detector, pixel_coordinates_in_detector_system
from .types import (
    ChopperDistance,
    ChopperFrequency,
    ChopperPhase,
    ChopperSeparation,
    RawChopper,
)


def load_detector(
    file_path: Filename[RunType], detector_name: NeXusDetectorName[RunType]
) -> LoadedNeXusDetector[RunType]:
    return nexus.load_detector(file_path=file_path, detector_name=detector_name)


def load_events(
    detector: LoadedNeXusDetector[RunType],
    detector_rotation: DetectorRotation[RunType],
    sample_rotation: SampleRotation[RunType],
    chopper_phase: ChopperPhase[RunType],
    chopper_frequency: ChopperFrequency[RunType],
    chopper_distance: ChopperDistance[RunType],
    chopper_separation: ChopperSeparation[RunType],
    sample_size: SampleSize[RunType],
) -> RawDetectorData[RunType]:
    detector_numbers = pixel_coordinates_in_detector_system()
    data = (
        nexus.extract_detector_data(detector)
        .bins.constituents["data"]
        .group(detector_numbers.data.flatten(to='event_id'))
        .fold("event_id", sizes=detector_numbers.sizes)
    )
    for name, coord in detector_numbers.coords.items():
        data.coords[name] = coord

    data.coords['z_index'] = (
        Detector.nWires * data.coords['blade'] + data.coords['wire']
    )

    if data.bins.constituents["data"].data.variances is None:
        data.bins.constituents["data"].data.variances = data.bins.constituents[
            "data"
        ].data.values

    data.coords["sample_rotation"] = sample_rotation.to(unit='rad')
    data.coords["detector_rotation"] = detector_rotation.to(unit='rad')
    data.coords["chopper_phase"] = chopper_phase
    data.coords["chopper_frequency"] = chopper_frequency
    data.coords["chopper_separation"] = sc.abs(chopper_separation)
    data.coords["L1"] = sc.abs(chopper_distance)
    data.coords["L2"] = data.coords['distance_in_detector'] + Detector.distance
    data.coords["sample_size"] = sample_size
    return RawDetectorData[RunType](data)


def amor_chopper(f: Filename[RunType]) -> RawChopper[RunType]:
    return next(load_nx(f, "NXentry/NXinstrument/NXdisk_chopper"))


def load_amor_chopper_distance(ch: RawChopper[RunType]) -> ChopperDistance[RunType]:
    # We know the value has unit 'mm'
    return sc.scalar(ch["distance"], unit="mm").to(unit="mm")


def load_amor_chopper_separation(ch: RawChopper[RunType]) -> ChopperSeparation[RunType]:
    # We know the value has unit 'mm'
    return sc.scalar(ch["pair_separation"], unit="mm").to(unit="mm")


def load_amor_ch_phase(ch: RawChopper[RunType]) -> ChopperPhase[RunType]:
    p = ch["phase"]["value"].coords["average_value"].value
    if getattr(p, "unit", None):
        return p
    raise ValueError("No unit was found for the chopper phase")


def load_amor_ch_frequency(ch: RawChopper[RunType]) -> ChopperFrequency[RunType]:
    f = ch["rotation_speed"]["value"].coords["average_value"]
    if getattr(f, "unit", None):
        return f
    raise ValueError("No unit was found for the chopper frequency")


def load_amor_sample_rotation(fp: Filename[RunType]) -> SampleRotation[RunType]:
    (mu,) = load_nx(fp, "NXentry/NXinstrument/master_parameters/mu")
    # Jochens Amor code reads the first value of this log
    # see https://github.com/jochenstahn/amor/blob/140e3192ddb7e7f28acee87e2acaee65ce1332aa/libeos/file_reader.py#L272  # noqa: E501
    # might have to change if this field ever becomes truly time-dependent
    return sc.scalar(mu['value'].data['dim_1', 0]['time', 0].value, unit='deg')


def load_amor_detector_rotation(fp: Filename[RunType]) -> DetectorRotation[RunType]:
    (nu,) = load_nx(fp, "NXentry/NXinstrument/master_parameters/nu")
    # Jochens Amor code reads the first value of this log
    # see https://github.com/jochenstahn/amor/blob/140e3192ddb7e7f28acee87e2acaee65ce1332aa/libeos/file_reader.py#L272  # noqa: E501
    # might have to change if this field ever becomes truly time-dependent
    return sc.scalar(nu['value'].data['dim_1', 0]['time', 0].value, unit='deg')


providers = (
    load_detector,
    load_events,
    load_amor_ch_frequency,
    load_amor_ch_phase,
    load_amor_chopper_distance,
    load_amor_chopper_separation,
    load_amor_sample_rotation,
    load_amor_detector_rotation,
    amor_chopper,
)
