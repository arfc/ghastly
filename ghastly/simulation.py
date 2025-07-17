import numpy as np

rng = np.random.default_rng()


class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''

    def __init__(self, r_pebble, t_final, pf,
                 recirc_target=1, recirc_rate=1, recirc_hz = 1,
                 v_center=1, v_mid=1, v_wall=1, 
                 core_intake={},core_main={}, core_outtake={}, recirc = {},
                 k_rate=0.001, down_flow=True, 
                 seed=rng.integers(1000000, 100000000)):
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
        core_intake : dict
            Dictionary containing key:value pairs where each key is the name
            of a core element within the core intake zone, and the value is a
            dictionary with key:value pairs corresponding to a core element
            parameter and its value.
        core_main : dict
            As core_intake, but for the main zone of the core.
        core_outtake : dict
            As core_intake, but for the outtake zone of the core.
        k_rate : float
            Contraction rate, to be used with OpenMC's pack_spheres function.
            Default is 0.001.
        down_flow : bool
            The direction of flow in the reactor core.  True means it is
            flowing downward, False means it is flowing upward.  Default is
            down, i.e., True.
        seed : int
            Seed used in LAMMPS simulations whenever a seed is required.
            Default is a random integer 6 to 7 digits long.
        '''
        self.r_pebble = r_pebble
        self.pebble_volume = (4/3)*np.pi*(r_pebble**3)
        self.t_final = t_final
        self.pf = pf
        self.recirc_target = recirc_target
        self.recirc_rate = recirc_rate
        self.recirc_hz = recirc_hz
        self.v_center = v_center
        self.v_mid = v_mid
        self.v_wall = v_wall
        self.core_intake = core_intake
        self.core_main = core_main
        self.core_outtake = core_outtake
        self.recirc = recirc
        self.k_rate = k_rate
        self.down_flow = down_flow
        self.seed = seed
