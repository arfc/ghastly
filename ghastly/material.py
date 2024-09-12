import numpy as np


class Mat():
    '''
    Class containing OpenMC material metadata necessary for pebble tracking
    '''

    def __init__(self, material, pass_num, reg_id, l_type=0):
        '''
        Initializes a single instance of the Mat class.

        Parameters
        ----------
        material : OpenMC Material
            The OpenMC material object this metadata applies to.
        pass_num : int
            Integer corresponding to the number of passes a pebble with this
            material would have completed.  A fresh pebble would have a
            pass_num of 0, for example.
        reg_id  : str
            reg_id for the region this material corresponds to.
        l_type : int
            The type in LAMMPS that corresponds to this material.  These will
            be determined by ghastly during simulation set up.
        '''
        self.material = material
        self.pass_num = pass_num
        self.reg_id = reg_id
        self.l_type = l_type
