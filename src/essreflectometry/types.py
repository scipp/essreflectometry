from typing import NewType, TypeVar

import sciline
import scipp as sc

Reference = NewType('Reference', str)
Sample = NewType('Sample', str)
Run = TypeVar('Run', Reference, Sample)


class RawData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data"""


class WavelengthData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Raw data transformed to wavelength"""


class ThetaData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Wavelength data transformed to theta"""


class QData(sciline.Scope[Run, sc.DataArray], sc.DataArray):
    """Theta data transformed to momentum transfer"""


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


NormalizedIofQ = NewType('NormalizedIofQ', sc.DataArray)
QResolution = NewType('QResolution', sc.Variable)


''' Parameters for the workflow '''

QBins = NewType('QBins', sc.Variable)


class Filename(sciline.Scope[Run, str], str):
    """Filename of the raw data"""


# TODO What do they mean?
# Supermirror parameters
MValue = NewType('MValue', sc.Variable)
CriticalEdge = NewType('CriticalEdge', sc.Variable)
Alpha = NewType('Alpha', sc.Variable)
SupermirrorParameter = TypeVar('SupermirrorParameter', MValue, CriticalEdge, Alpha)


class SampleRotation(sciline.Scope[Run, sc.Variable], sc.Variable):
    """The rotation of the sample / the reference sample."""


class BeamSize(sciline.Scope[Run, sc.Variable], sc.Variable):
    """FWHM of the neutron beam."""


class DetectorSpatialResolution(sciline.Scope[Run, sc.Variable], sc.Variable):
    # TODO
    """Don't know what this is."""


class SampleSize(sciline.Scope[Run, sc.Variable], sc.Variable):
    # TODO is this radius or total length?
    """Size of the sample."""


Gravity = NewType('Gravity', sc.Variable)
