import numpy as np
import pytest
from ghastly import pebble


def test_Pebble():
    '''
    Tests the pebble class.
    '''
    test_peb = pebble.Pebble(10, np.array([1.0, 2.0, 3.0]), 13,
                             1, 4, 3, 0)

    assert list(test_peb.coords) == [1.0, 2.0, 3.0]
    assert test_peb.uid == 10
    assert test_peb.recirc == 0
    assert test_peb.zone == 1
    assert test_peb.velocity == 13
