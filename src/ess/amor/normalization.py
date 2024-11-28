import scipp as sc

from ..reflectometry.conversions import reflectometry_q
from ..reflectometry.supermirror import (
    Alpha,
    CriticalEdge,
    MValue,
    supermirror_reflectivity,
)
from ..reflectometry.types import (
    QBins,
    QResolution,
    ReducedReference,
    ReducibleData,
    Reference,
    ReferenceRun,
    ReflectivityOverQ,
    ReflectivityOverZW,
    SampleRun,
    WavelengthBins,
    ZIndexLimits,
)
from .conversions import angle_of_reflection


def _add_pre_reduction_masks(da, zindex_limits):
    da.masks['z_range'] = (da.coords["iz"] < zindex_limits[0]) | (
        da.coords["iz"] > zindex_limits[1]
    )


def reduce_reference(
    reference: ReducibleData[ReferenceRun],
    wavelength_bins: WavelengthBins,
    critical_edge: CriticalEdge,
    mvalue: MValue,
    alpha: Alpha,
) -> ReducedReference:
    R = supermirror_reflectivity(
        reflectometry_q(
            wavelength=reference.bins.coords['wavelength'],
            angle_of_reflection=reference.bins.coords['angle_of_reflection'],
        ),
        c=critical_edge,
        m=mvalue,
        alpha=alpha,
    )
    reference.bins.masks['invalid'] = sc.isnan(R)
    reference /= R
    return reference.bins.concat(('stripe',)).hist(wavelength=wavelength_bins)


def reduce_sample_over_q(
    sample: ReducibleData[SampleRun],
    reference: Reference,
    qbins: QBins,
    wavelength_bins: WavelengthBins,
    zlims: ZIndexLimits,
    qresolution: QResolution,
) -> ReflectivityOverQ:
    sample.bins.coords['Q'] = reflectometry_q(
        sample.bins.coords['wavelength'],
        sample.bins.coords['angle_of_reflection'],
    )
    sample.bins.masks['wavelength'] = (
        wavelength_bins[0] > sample.bins.coords['wavelength']
    ) | (wavelength_bins[-1] < sample.bins.coords['wavelength'])
    _add_pre_reduction_masks(sample, zlims)
    R = sample.bins.concat().bin(Q=qbins) / reference.flatten(to='Q').hist(Q=qbins).data
    R.coords['Q_resolution'] = qresolution.data
    return R


def reduce_sample_over_zw(
    sample: ReducibleData[SampleRun],
    reference: Reference,
    wavelength_bins: WavelengthBins,
    zlims: ZIndexLimits,
) -> ReflectivityOverZW:
    sample.bins.coords['Q'] = reflectometry_q(
        sample.bins.coords['wavelength'],
        sample.bins.coords['angle_of_reflection'],
    )
    sample.bins.masks['wavelength'] = (
        wavelength_bins[0] > sample.bins.coords['wavelength']
    ) | (wavelength_bins[-1] < sample.bins.coords['wavelength'])
    _add_pre_reduction_masks(sample, zlims)
    R = sample.bins.concat(('stripe',)).bin(wavelength=wavelength_bins) / reference.data
    R.masks["too_few_events"] = reference.data < sc.scalar(1, unit="counts")
    return R


def evaluate_reference(
    reference: ReducedReference,
    sample: ReducibleData[SampleRun],
    qbins: QBins,
    zlims: ZIndexLimits,
) -> Reference:
    ref = reference.copy()
    ref.coords['angle_of_reflection'] = angle_of_reflection(
        wavelength=sc.midpoints(ref.coords["wavelength"]),
        divergence_angle=ref.coords["pixel_divergence_angle"],
        L2=ref.coords["L2"],
        sample_rotation=sample.coords["sample_rotation"],
        detector_rotation=sample.coords["detector_rotation"],
    )
    ref.coords['Q'] = reflectometry_q(
        wavelength=sc.midpoints(ref.coords["wavelength"]),
        angle_of_reflection=ref.coords["angle_of_reflection"],
    )
    _add_pre_reduction_masks(ref, zlims)
    return sc.values(ref)


providers = (
    reduce_reference,
    reduce_sample_over_q,
    reduce_sample_over_zw,
    evaluate_reference,
)
