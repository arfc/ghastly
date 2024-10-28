import pytest
from ghastly import read_input
from ghastly import simulation

def test_InputBlock():
    '''
    Test InputBlock class
    '''
    test_input = read_input.InputBlock("sample_input.json")

    assert test_input.sim_var["r_pebble"] == 0.03
    assert test_input.core_intake_var["intake_cone"]["z_max"] == 1.1
    assert test_input.core_main_var["main_cyl"]["type"] == "cylinder"
    assert test_input.core_outtake_var["outtake_cyl"]["r"] == 0.1
    assert test_input.lammps_var["density_rho"] == 1700

def test_create_obj():
    '''
    Test InputBlock's create_obj method
    '''
    test_input = read_input.InputBlock("sample_input.json")
    test_simblock = test_input.create_obj()

    assert type(test_simblock) == simulation.Sim

def test_create_core_zone():
    '''
    Test InputBlock's create_core_zone method
    '''
    pass

def create_sim_block():
    '''
    Tests InputBlock's create_sim_block method
    '''
    pass
