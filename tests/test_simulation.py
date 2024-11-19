import pytest
import numpy as np
from ghastly import simulation


def test_Sim():
    '''
    Test the Sim class.
    '''
    test_sim = simulation.Sim(1.0, 10.0, 0.6)
    assert test_sim.r_pebble == 1.0
    assert test_sim.pebble_volume == pytest.approx((4/3)*np.pi*1.0**3)
    assert test_sim.t_final == 10.0
    assert test_sim.pf == 0.6
    assert test_sim.core_intake == {}
    assert test_sim.core_main == {}
    assert test_sim.core_outtake == {}
    assert test_sim.k_rate == 0.001
    assert test_sim.down_flow == True
    assert test_sim.seed >= 1000000 or test_sim.seed <= 100000000
