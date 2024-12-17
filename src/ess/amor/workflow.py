from ..reflectometry.corrections import correct_by_footprint, correct_by_proton_current
from ..reflectometry.types import (
    BeamDivergenceLimits,
    ProtonCurrent,
    RawDetectorData,
    ReducibleData,
    RunType,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)
from .conversions import add_coords, add_masks


def add_coords_masks_and_apply_corrections(
    da: RawDetectorData[RunType],
    ylim: YIndexLimits,
    zlims: ZIndexLimits,
    bdlim: BeamDivergenceLimits,
    wbins: WavelengthBins,
    proton_current: ProtonCurrent[RunType],
) -> ReducibleData[RunType]:
    """
    Computes coordinates, masks and corrections that are
    the same for the sample measurement and the reference measurement.
    """
    da = add_coords(da, proton_current)
    da = add_masks(da, ylim, zlims, bdlim, wbins)
    correct_by_footprint(da)
    if 'proton_current' in da.bins.coords:
        correct_by_proton_current(da)
    return da


providers = (add_coords_masks_and_apply_corrections,)
