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


QResolution = NewType('QResolution', sc.Variable)


class FootprintCorrected(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Experiment data corrected by footprint on sample"""


SpecularReflectionCoordTransformGraph = NewType(
    'SpecularReflectionCoordTransformGraph', dict
)

CalibratedReference = NewType('CalibratedReference', sc.DataArray)


class Histogrammed(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Histogrammmed by Q and detector_number"""


class Normalized(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Normalized histogram"""


NormalizedIOfQ = NewType('NormalizedIOfQ', sc.DataArray)


''' Parameters for the workflow '''

QBins = NewType('QBins', sc.Variable)


class Filename(sciline.Scope[Run, str], str):
    """Filename of the raw data"""
