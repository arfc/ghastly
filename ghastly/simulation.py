import numpy as np
from jinja2 import Environment, FileSystemLoader

rng = np.random.default_rng()

class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''
    def __init__(self, r_pebble, t_final, pf, core_intake = {},
                 core_main = {}, core_outtake = {}, k_rate = 0.001, 
                 down_flow=True, seed = rng.integers(1000000,100000000)):
        '''
        Initializes the Sim class.

        Parameters
        ----------
        r_pebble : float
            Radius of the pebbles in the simulation [m].
        t_final : float
            Total reactor-time that is being simulated [s].
        pf : float
            Target packing fraction in core region.  Determines total number
            of pebbles, but packing fraction may not match this value in
            all areas due to settling.
        down_flow : bool
            The direction of flow in the reactor core.  True means it is
            flowing downward, False means it is flowing upward.  Default is
            down.
        '''
        self.r_pebble = r_pebble
        self.pebble_volume = (4/3)*np.pi*(r_pebble**3)
        self.t_final = t_final
        self.pf = pf
        self.core_intake = core_intake
        self.core_main = core_main
        self.core_outtake = core_outtake
        self.k_rate = k_rate
        self.down_flow = down_flow
        self.seed = seed


                   
