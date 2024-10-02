from typing import NewType, TypeVar

import sciline
import scipp as sc

ReferenceRun = NewType("ReferenceRun", str)
SampleRun = NewType("SampleRun", str)
RunType = TypeVar("RunType", ReferenceRun, SampleRun)


class NeXusDetectorName(sciline.Scope[RunType, str], str):
    """Name of the detector in the nexus file containing the events of the RunType"""


class DetectorPosition(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """Positions of the detector pixels, relative to the source(?), as a 3d-vector"""


class SamplePosition(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """The position of the sample relative to the source(?)."""


class IncidentBeam(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """Incident beam vector."""


class SpecularReflectionCoordTransformGraph(sciline.Scope[RunType, dict], dict):
    """Coordinate transformation graph for specular reflection"""


class RawDetectorData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Event time data from nexus file,
    binned by `detector_number` (pixel of the detector frame)."""


class EventData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Event data from nexus file, binned by pulse.
    Holds time dependent masks and coords."""


class RawEventData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Raw event data from nexus file, binned by pulse."""


class LoadedNeXusDetector(sciline.Scope[RunType, sc.DataGroup], sc.DataGroup):
    """NXdetector loaded from file"""


class ReducibleDetectorData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Event time data after correcting tof, ready for reduction"""


class DataWithScatteringCoordinates(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Event data with added coordinates such as incident angle (theta),
    wavelength, and momentum transfer (Q)"""


class MaskedData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Event data that has been masked in wavelength and logical detector coordinates"""


class FootprintCorrection(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Weights of the footprint of the beam on the sample."""


class ProtonCurrentCorrection(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """Weights of the proton current correction"""


class CorrectedData(sciline.Scope[RunType, sc.DataArray], sc.DataArray):
    """MaskedData corrected by footprint and proton current"""


ReferenceIntensity = NewType("ReferenceIntensity", sc.DataArray)
"""Intensity distribution of the reference measurement in (z, wavelength)"""

IdealReferenceIntensity = NewType("IdealReferenceIntensity", sc.DataArray)
"""Intensity distribution on the detector for a sample with :math`R(Q) = 1`"""

NormalizationFactor = NewType("NormalizationFactor", sc.DataArray)
""":code`IdealReferenceIntensity` with added coordinate "sample"-Q"""

NormalizedIofQ = NewType("NormalizedIofQ", sc.DataArray)
"""Intensity histogram over momentum transfer
normalized by the calibrated reference measurement."""

ReflectivityData = NewType("ReflectivityData", sc.DataArray)
"""Reflectivity "per event". Event data weighted by the expected
intensity at the coordinates of the event."""

QResolution = NewType("QResolution", sc.Variable)
"""Resolution term for the momentum transfer for each bin of QBins."""


""" Parameters for the workflow """

QBins = NewType("QBins", sc.Variable)
"""Bins for the momentum transfer histogram."""

WavelengthBins = NewType("WavelengthBins", sc.Variable)
"""Bins for the wavelength histogram, also used to filter the event data."""


class Filename(sciline.Scope[RunType, str], str):
    """Filename of an event data nexus file."""


class SampleRotation(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """The rotation of the sample relative to the center of the incoming beam."""


class DetectorRotation(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """The rotation of the detector relative to the horizon"""


class BeamSize(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """Full-Width-Half-maximum of the incoming beam."""


class DetectorSpatialResolution(sciline.Scope[RunType, sc.Variable], sc.Variable):
    # TODO what is the definition of this?
    """Spatial resolution of the detector."""


class SampleSize(sciline.Scope[RunType, sc.Variable], sc.Variable):
    # TODO is this radius or total length?
    """Size of the sample. If None it is assumed to be the same as the reference."""


class ProtonCurrent(sciline.Scope[RunType, sc.Variable], sc.Variable):
    """Proton current log of the measurement"""


Gravity = NewType("Gravity", sc.Variable)
"""This parameter determines if gravity is taken into account
when computing the scattering angle and momentum transfer."""


YIndexLimits = NewType("YIndexLimits", tuple[sc.Variable, sc.Variable])
"""Limit of the (logical) 'y' detector pixel index"""


ZIndexLimits = NewType("ZIndexLimits", tuple[sc.Variable, sc.Variable])
"""Limit of the (logical) 'z' detector pixel index"""


BeamDivergenceLimits = NewType("BeamDivergenceLimits", tuple[sc.Variable, sc.Variable])
"""Limit of the beam divergence"""


ReferenceFilePath = NewType("ReferenceFilePath", str)
"""Path to the cached normalization matrix"""
