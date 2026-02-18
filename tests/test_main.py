import pytest
import numpy as np
import glob
import os
from ghastly import main
from ghastly import read_input

test_input = read_input.InputBlock("fill_input.json")
test_sim = test_input.create_obj()

def test_pack_core():
    '''
    Tests the main function pack_cyl.
    '''
    test_cyl = test_sim.core_main["main_cyl"]
    test_coords = main._pack_core(test_sim, 
                                  test_cyl,
                                  rough_pf=0.2, 
                                  crp = False, 
                                  openmc = False)
    stack = np.vstack(test_coords)
    rlim = test_cyl.r - test_sim.r_pebble
    zlim_min = test_cyl.zmin + test_sim.r_pebble
    zlim_max = test_cyl.zmax - test_sim.r_pebble
    assert stack.min(axis=0)[0] > test_cyl.x_c - rlim
    assert stack.min(axis=0)[1] > test_cyl.y_c - rlim
    assert stack.min(axis=0)[2] > zlim_min
    assert stack.max(axis=0)[0] < test_cyl.x_c + rlim
    assert stack.max(axis=0)[1] < test_cyl.y_c + rlim
    assert stack.max(axis=0)[2] < zlim_max


def test_find_box_bounds():
    '''
    Tests the find_box_bounds function in main.py.
    '''
    bound_limits = main.find_box_bounds(test_sim)

    assert bound_limits['xb_min'] == pytest.approx(0.05 -0.5 - 5*0.03)
    assert bound_limits['yb_max'] == pytest.approx(0.5 + 5*0.03)
    assert bound_limits['zb_min'] == pytest.approx(-0.1 - 5*0.03)
    assert bound_limits['zb_max'] == pytest.approx(1.2 + 5*0.03)


def test_write_variable_block():
    '''
    Tests _write_variable_block in main.py.
    '''
    variables = main._write_variable_block(test_input, test_sim,
                                           "variables.txt", 
                                           "lammps/variable_template.txt")

    assert variables['r_vessel'] == pytest.approx(0.5)
    assert variables['vessel_zmin'] == pytest.approx(-0.1)
    assert variables['vessel_zmax'] == pytest.approx(1.0)
    assert variables['vessel_x_c'] == pytest.approx(0.05)
    assert variables['vessel_y_c'] == pytest.approx(0.0)

def test_fill_core():
    '''
    Tests fill_core in main.py.
    '''
    main.fill_core("fill_input.json", rough_pf = 0.15)

def test_recirc_pebbles():
    '''
    Tests recirc_pebbles in main.py
    '''
    main.recirc_pebbles("f1_input.json", "rough-pack.txt", 
                        "f1_main.txt", "lammps/f1main_template.txt")
    main.recirc_pebbles("f2_input.json", "rough-pack.txt", 
                        "f2_main.txt", "lammps/f2main_template.txt")


def test_cleanup():
    '''
    Removes files after testing
    '''
    for f in glob.glob("*.txt"):
        os.remove(f)





