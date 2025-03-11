# Copyright (c) 2025 Scipp contributors (https://github.com/scipp)
from ..reflectometry.conversions import (
    add_coords,
    add_proton_current_coord,
    add_proton_current_mask,
)
from ..reflectometry.corrections import correct_by_proton_current
from ..reflectometry.types import (
    BeamDivergenceLimits,
    CoordTransformationGraph,
    ProtonCurrent,
    RawDetectorData,
    ReducibleData,
    RunType,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)
from .corrections import correct_by_footprint
from .maskings import add_masks


def add_coords_masks_and_apply_corrections(
    da: RawDetectorData[RunType],
    ylim: YIndexLimits,
    zlims: ZIndexLimits,
    bdlim: BeamDivergenceLimits,
    wbins: WavelengthBins,
    proton_current: ProtonCurrent[RunType],
    graph: CoordTransformationGraph,
) -> ReducibleData[RunType]:
    """
    Computes coordinates, masks and corrections that are
    the same for the sample measurement and the reference measurement.
    """
    da = add_coords(da, graph)
    da = add_masks(da, ylim, zlims, bdlim, wbins)
    da = correct_by_footprint(da)

    if len(proton_current) != 0:
        da = add_proton_current_coord(da, proton_current)
        da = add_proton_current_mask(da)
        da = correct_by_proton_current(da)

    return ReducibleData[RunType](da)


providers = (add_coords_masks_and_apply_corrections,)
