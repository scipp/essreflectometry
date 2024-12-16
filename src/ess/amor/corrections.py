import scipp as sc

from ..reflectometry.corrections import footprint_on_sample


def correct_by_footprint(da: sc.DataArray) -> None:
    da /= footprint_on_sample(
        da.bins.coords['theta'],
        da.coords['beam_size'],
        da.coords['sample_size'],
    )
