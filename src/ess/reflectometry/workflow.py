from collections.abc import Hashable, Sequence
from itertools import chain

import pandas as pd
import sciline
import scipp as sc

from ess.amor.types import RawChopper
from ess.reflectometry.orso import (
    OrsoExperiment,
    OrsoOwner,
    OrsoSample,
    OrsoSampleFilenames,
)
from ess.reflectometry.types import (
    Filename,
    FootprintCorrectedData,
    RunType,
    SampleRotation,
    SampleRun,
)


def _set_first_if_lost(coord, da, das):
    if coord not in da.coords:
        da.coords[coord] = das[0].coords[coord]


def _concatenate_event_lists(*das):
    da = sc.reduce(das).bins.concat()
    _set_first_if_lost('position', da, das)
    _set_first_if_lost('sample_rotation', da, das)
    _set_first_if_lost('detector_rotation', da, das)
    return da


def with_filenames(
    workflow, runtype: Hashable, runs: Sequence[Filename[RunType]]
) -> sciline.Pipeline:
    axis_name = f'{str(runtype).lower()}_runs'
    df = pd.DataFrame({Filename[runtype]: runs}).rename_axis(axis_name)
    wf = workflow.copy()

    mapped = wf.map(df)

    wf[FootprintCorrectedData[runtype]] = mapped[
        FootprintCorrectedData[runtype]
    ].reduce(index=axis_name, func=_concatenate_event_lists)
    wf[RawChopper[runtype]] = mapped[RawChopper[runtype]].reduce(
        index=axis_name, func=lambda x, *_: x
    )
    wf[SampleRotation[runtype]] = mapped[SampleRotation[runtype]].reduce(
        index=axis_name, func=lambda x, *_: x
    )

    if runtype is SampleRun:
        wf[OrsoSample] = mapped[OrsoSample].reduce(
            index=axis_name, func=lambda x, *_: x
        )
        wf[OrsoExperiment] = mapped[OrsoExperiment].reduce(
            index=axis_name, func=lambda x, *_: x
        )
        wf[OrsoOwner] = mapped[OrsoOwner].reduce(index=axis_name, func=lambda x, *_: x)
        wf[OrsoSampleFilenames] = mapped[OrsoSampleFilenames].reduce(
            index=axis_name, func=lambda *x: list(chain(*x))
        )
    return wf
