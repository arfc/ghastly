import numpy as np
from ghastly import core
import pytest


def test_CylCore():
    '''
    Tests the CylCore class
    '''

    test_cyl = core.CylCore(0.5, 1.0, 2.0, 4.0, 3.0)
    assert test_cyl.r == 0.5
    assert test_cyl.x_c == 1.0
    assert test_cyl.y_c == 2.0
    assert test_cyl.z_max == 4.0
    assert test_cyl.z_min == 3.0
    assert test_cyl.h == 1.0
    assert test_cyl.volume == pytest.approx(np.pi * (0.5**2) * 1.0)


def test_ConeCore():
    '''
    Tests the ConeCone class
    '''
    test_cone = core.ConeCore(1.5, 0.5, 1.0, 2.0, 4.0, 3.0)
    assert test_cone.r_upper == 1.5
    assert test_cone.r_lower == 0.5
    assert test_cone.x_c == 1.0
    assert test_cone.y_c == 2.0
    assert test_cone.z_max == 4.0
    assert test_cone.z_min == 3.0
    assert test_cone.h == 1.0
    assert test_cone.volume == pytest.approx((1 / 3) * np.pi * 1.0
                                             * (1.5**2 + 0.5**2 + 1.5 * 0.5))
