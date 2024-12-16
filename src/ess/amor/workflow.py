from ..reflectometry.types import (
    BeamDivergenceLimits,
    RawDetectorData,
    ReducibleData,
    RunType,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)
from .conversions import add_coords, add_masks
from .corrections import correct_by_footprint


def add_coords_masks_and_apply_corrections(
    da: RawDetectorData[RunType],
    ylim: YIndexLimits,
    zlims: ZIndexLimits,
    bdlim: BeamDivergenceLimits,
    wbins: WavelengthBins,
) -> ReducibleData[RunType]:
    da = add_coords(da)
    da = add_masks(da, ylim, zlims, bdlim, wbins)
    correct_by_footprint(da)
    return da


providers = (add_coords_masks_and_apply_corrections,)
