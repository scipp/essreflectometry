import scipp as sc

from .tools import fwhm_to_std


def footprint_on_sample(
    theta: sc.Variable,
    beam_size: sc.Variable,
    sample_size: sc.Variable,
) -> sc.Variable:
    """
    Computes the fraction of the beam hitting the sample.
    Depends on :math:`\\theta`.

    Parameters
    ----------
    theta:
        Incidence angle relative to sample surface.
    beam_size:
        Full width half maximum of the beam.
    sample_size:
        Size of the sample.
        TODO: check what sample size actually means. Is it the sample diameter? etc.

    Returns
    -------
    :
       Fraction of beam hitting the sample.
    """
    size_of_beam_on_sample = beam_size / sc.sin(theta)
    return sc.erf(fwhm_to_std(sample_size / size_of_beam_on_sample))
