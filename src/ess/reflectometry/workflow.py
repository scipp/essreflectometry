from collections.abc import Hashable, Sequence

import pandas as pd
import sciline
import scipp as sc

from ess.reflectometry.types import Filename, FootprintCorrectedData, RunType


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

    wf[FootprintCorrectedData[runtype]] = (
        wf[FootprintCorrectedData[runtype]]
        .map(df)
        .reduce(index=axis_name, func=_concatenate_event_lists)
    )
    return wf
