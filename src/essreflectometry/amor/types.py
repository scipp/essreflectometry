from typing import NewType, TypeVar

import sciline
import scipp as sc

from ..types import Run

WavelengthResolution = NewType('WavelengthResolution', sc.Variable)
AngularResolution = NewType('AngularResolution', sc.Variable)
SampleSizeResolution = NewType('SampleSizeResolution', sc.Variable)


class BeamlineParams(sciline.Scope[Run, dict], dict):
    """Parameters describing the beamline"""


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


class ChopperFrequency(sciline.Scope[Run, sc.Variable], sc.Variable):
    """Frequency of the choppers used in the run."""


class ChopperPhase(sciline.Scope[Run, sc.Variable], sc.Variable):
    """Phase of the choppers in the run."""


class Chopper1Position(sciline.Scope[Run, sc.Variable], sc.Variable):
    """Position of the first chopper relative the source of the beam."""


class Chopper2Position(sciline.Scope[Run, sc.Variable], sc.Variable):
    """Position of the second chopper relative to the source of the beam."""
