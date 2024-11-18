import pytest
import numpy as np
from ghastly import ghastly
from ghastly import read_input

def test_pack_cyl():
    '''
    Tests the ghastly function pack_cyl
    '''
    test_input = read_input.InputBlock("sample_input.json")
    test_sim = test_input.create_obj()
    test_cyl = test_sim.core_main["main_cyl"]
    test_coords = ghastly.pack_cyl(test_sim, test_sim.core_main["main_cyl"],
                                   rough_pf = 0.2)
    n_pebs = int((0.2*test_cyl.volume)/test_sim.pebble_volume)
    assert len(test_coords) == pytest.approx(n_pebs, abs=1.0)
    stack = np.vstack(test_coords)
    assert stack.min(axis=0)[0] > test_cyl.x_c - test_cyl.r
    assert stack.min(axis=0)[1] > test_cyl.y_c - test_cyl.r
    assert stack.min(axis=0)[2] > test_cyl.z_min
    assert stack.max(axis=0)[0] < test_cyl.x_c + test_cyl.r
    assert stack.max(axis=0)[1] < test_cyl.y_c + test_cyl.r
    assert stack.max(axis=0)[2] < test_cyl.z_max

def test_find_box_bounds():
    '''
    Tests the find_box_bounds function in ghastly.py
    '''
    test_input = read_input.InputBlock("sample_input.json")
    test_sim = test_input.create_obj()
    x_b, y_b, z_b = ghastly.find_box_bounds(test_sim)

    assert x_b["low"] == pytest.approx(0.05 - 1.05*0.5)
    assert y_b["up"] == pytest.approx(0.0 + 1.05*0.5)
    assert z_b["low"] == pytest.approx(0.55 - 1.05*.65)


