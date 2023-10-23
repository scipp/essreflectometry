from typing import NewType, TypeVar

import sciline
import scipp as sc

Reference = NewType('Reference', str)
Sample = NewType('Sample', str)
Run = TypeVar('Run', Reference, Sample)


class Raw(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data"""


class WavelengthData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data transformed to wavelength"""


class ThetaData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Wavelength data transformed to theta"""


class QData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Theta data transformed to momentum transfer"""


QStd = NewType('QStd', sc.Variable)


class FootprintCorrected(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Experiment data corrected by footprint on sample"""


WavelengthResolution = NewType('WavelengthResolution', sc.Variable)
AngularResolution = NewType('AngularResolution', sc.Variable)
SampleSizeResolution = NewType('SampleSizeResolution', sc.Variable)


class BeamlineParams(sciline.Scope[Run, dict], dict):
    """Parameters describing the beamline"""


SpecularReflectionCoordTransformGraph = NewType(
    'SpecularReflectionCoordTransformGraph', dict
)

CalibratedReference = NewType('CalibratedReference', sc.DataArray)
HistogramContent = TypeVar(
    'HistogramContent', Sample, Reference, CalibratedReference, AngularResolution
)


class Histogrammed(sciline.Scope[HistogramContent, sc.DataArray], sc.DataArray):
    """Histogrammmed by Q and detector_number"""


CalibratedRun = TypeVar('CalibratedRun', CalibratedReference, Sample)


class Normalized(sciline.Scope[CalibratedRun, sc.DataArray], sc.DataArray):
    """Normalized histogram"""


NormalizedIOverQ = NewType('NormalizedIOverQ', sc.DataArray)


''' Parameters for the workflow '''

QBins = NewType('QBins', sc.Variable)


class SampleRotation(sciline.Scope[Run, sc.Variable], sc.Variable):
    """The rotation of the sample / the reference sample"""


class Filename(sciline.Scope[Run, str], str):
    """Filename of the raw data"""
