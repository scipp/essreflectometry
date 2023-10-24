from typing import NewType, TypeVar

import sciline
import scipp as sc

from ..reflectometry.types import Run

MValue = NewType('MValue', str)
CriticalEdge = NewType('CriticalEdge', str)
Alpha = NewType('Alpha', str)
SupermirrorParameter = TypeVar('SupermirrorParameter', MValue, CriticalEdge, Alpha)


class Supermirror(sciline.Scope[SupermirrorParameter, sc.Variable], sc.Variable):
    """Supermirror parameter scope"""


class SampleRotation(sciline.Scope[Run, sc.Variable], sc.Variable):
    """The rotation of the sample / the reference sample"""


class BeamSize(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


class DetectorSpatialResolution(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


class SampleSize(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


Gravity = NewType('Gravity', sc.Variable)


class ChopperFrequency(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


class ChopperPhase(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


class Chopper1Position(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""


class Chopper2Position(sciline.Scope[Run, sc.Variable], sc.Variable):
    """parameter scope"""
