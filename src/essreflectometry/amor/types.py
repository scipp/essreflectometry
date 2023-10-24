from typing import NewType, TypeVar

import sciline
import scipp as sc

MValue = NewType('MValue', str)
CriticalEdge = NewType('CriticalEdge', str)
Alpha = NewType('Alpha', str)
SupermirrorParameter = TypeVar('SupermirrorParameter', MValue, CriticalEdge, Alpha)


class Supermirror(sciline.Scope[SupermirrorParameter, sc.Variable], sc.Variable):
    """Supermirror parameter scope"""
