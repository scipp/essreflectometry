import numpy as np
import pytest
import scipp as sc
from scipp.testing import assert_allclose

from ess.estia.calibration import solve_for_calibration_parameters


def generate_valid_calibration_parameters():
    I0 = np.random.random()
    Pp = np.random.random()
    Pa = -np.random.random()
    Ap = np.random.random()
    Aa = -np.random.random()
    Rspp = np.random.random()
    Rsaa = Rspp * np.random.random()
    return tuple(map(sc.scalar, (I0, Pp, Pa, Ap, Aa, Rspp, Rsaa)))


def intensity_from_parameters(I0, Pp, Pa, Ap, Aa, Rpp, Rpa, Rap, Raa):
    return (
        I0
        * (
            Rpp * (1 + Ap) * (1 + Pp)
            + Rpa * (1 - Ap) * (1 + Pp)
            + Rap * (1 + Ap) * (1 - Pp)
            + Raa * (1 - Ap) * (1 - Pp)
        ),
        I0
        * (
            Rpp * (1 + Aa) * (1 + Pp)
            + Rpa * (1 - Aa) * (1 + Pp)
            + Rap * (1 + Aa) * (1 - Pp)
            + Raa * (1 - Aa) * (1 - Pp)
        ),
        I0
        * (
            Rpp * (1 + Ap) * (1 + Pa)
            + Rpa * (1 - Ap) * (1 + Pa)
            + Rap * (1 + Ap) * (1 - Pa)
            + Raa * (1 - Ap) * (1 - Pa)
        ),
        I0
        * (
            Rpp * (1 + Aa) * (1 + Pa)
            + Rpa * (1 - Aa) * (1 + Pa)
            + Rap * (1 + Aa) * (1 - Pa)
            + Raa * (1 - Aa) * (1 - Pa)
        ),
    )


@pytest.mark.parametrize("seed", range(10))
def test_calibration_solve_recovers_input(seed):
    np.random.seed(seed)
    I0, Pp, Pa, Ap, Aa, Rspp, Rsaa = generate_valid_calibration_parameters()
    Io = intensity_from_parameters(
        I0, Pp, Pa, Ap, Aa, sc.scalar(1), sc.scalar(0), sc.scalar(0), sc.scalar(1)
    )
    Is = intensity_from_parameters(
        I0, Pp, Pa, Ap, Aa, Rspp, sc.scalar(0), sc.scalar(0), Rsaa
    )
    tuple(
        map(
            assert_allclose,
            solve_for_calibration_parameters(Io, Is),
            (I0, Pp, Pa, Ap, Aa, Rspp, Rsaa),
        )
    )
