from typing import NewType, TypeVar

import sciline
import scipp as sc

Reference = NewType('Reference', str)
Sample = NewType('Sample', str)
Run = TypeVar('Run', Reference, Sample)


class Filename(sciline.Scope[Run, str], str):
    """Filename of the raw data"""


class Raw(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data"""


class WavelengthData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data transformed to wavelength"""


class ThetaData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Wavelength data transformed to theta"""


class QData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Theta data transformed to momentum transfer"""


class FootprintCorrected(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Experiment data corrected by footprint on sample"""


CalibratedReference = NewType('CalibratedReference', sc.DataArray)
Resolutions = NewType('Resolutions', dict)
Normalized = NewType('Normalized', sc.DataArray)
CountsByMomentumTransfer = NewType('CountsByMomentumTransfer', sc.DataArray)

''' Parameters for the workflow '''
QBins = NewType('QBins', sc.Variable)
WavelengthBins = NewType('WavelengthBins', sc.Variable)
ThetaBins = NewType('ThetaBins', sc.Variable)


class SampleRotation(sciline.Scope[Run, sc.Variable], sc.Variable):
    """The rotation of the sample / the reference sample"""


class BeamlineParams(sciline.Scope[Run, dict], dict):
    """Parameters describing the beamline"""


SpecularReflectionCoordTransformGraph = NewType(
    'SpecularReflectionCoordTransformGraph', dict
)

QDataWithResolutions = NewType('QDataWithResolutions', sc.DataArray)
CorrectedQData = TypeVar('CorrectedQData', QData[Reference], QDataWithResolutions)


class HistogrammedByQ(sciline.Scope[CorrectedQData, sc.DataArray], sc.DataArray):
    """Histogrammmed by Q. Either reference data or sample data with resolutions."""


Normalizable = TypeVar(
    'Normalizable', HistogrammedByQ[QDataWithResolutions], CalibratedReference
)


class NormalizedData(sciline.Scope[Normalizable, sc.DataArray], sc.DataArray):
    """Normalized histogramm by Q."""


NormalizedIOverQ = NewType('NormalizedIOverQ', sc.DataArray)
