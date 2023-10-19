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


class Experiment(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Experiment data with added coordinates:
    wavelength, incidence angle, and momentum transfer"""


class FootprintCorrected(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Experiment data corrected by footprint on sample"""


CalibratedReference = NewType('CalibratedReference', sc.DataArray)
WithQResolution = NewType('WithQResolution', sc.DataArray)
Normalized = NewType('Normalized', sc.DataArray)
CountsByMomentumTransfer = NewType('CountsByMomentumTransfer', sc.DataArray)

''' Parameters for the workflow '''
MomentumTransferBins = NewType('MomentumTransferBins', sc.Variable)
WavelengthBins = NewType('WavelengthBins', sc.Variable)
ThetaBins = NewType('ThetaBins', sc.Variable)


class Rotation(sciline.Scope[Run, sc.Variable], sc.Variable):
    """The rotation of the sample / the reference sample"""


SpecularReflectionCoordTransformGraph = NewType(
    'SpecularReflectionCoordTransformGraph', dict
)
