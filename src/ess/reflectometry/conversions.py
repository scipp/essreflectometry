# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import scipp as sc
from scipp.constants import pi
from scippneutron._utils import elem_dtype


def reflectometry_q(
    wavelength: sc.Variable, angle_of_reflection: sc.Variable
) -> sc.Variable:
    """
    Compute momentum transfer from reflection angle.

    Parameters
    ----------
    wavelength:
        Wavelength values for the events.
    angle_of_reflection:
        Angle of reflection for the events.

    Returns
    -------
    :
        Q-values.
    """
    dtype = elem_dtype(wavelength)
    c = (4 * pi).astype(dtype)
    return c * sc.sin(angle_of_reflection.astype(dtype, copy=False)) / wavelength


providers = ()
