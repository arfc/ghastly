import numpy as np
import pytest
from ghastly import pebble

def test_Pebble():
    '''
    Tests the pebble class.
    '''
    test_peb = pebble.Pebble(np.array([1.0,2.0,3.0]), "test_reg", 
                                     0, 1, 27)

    assert list(test_peb.coords) == [1.0,2.0,3.0]
    assert test_peb.reg_id == "test_reg"
    assert test_peb.pass_num == 0
    assert test_peb.l_type == 1
    assert test_peb.pebble_id == 27

