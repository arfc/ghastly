import numpy as np

rng = np.random.default_rng()


class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''

    def __init__(self, r_pebble, t_final, pf, fidelity = 1,
                 recirc_target=1, recirc_hz = 1, 
                 core_intake={},core_main={}, core_outtake={}, recirc = {},
                 k_rate=0.001, down_flow=True, 
                 seed=rng.integers(1000000, 9999999)):
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
        fidelity : int
            Integer from the following list: [1, 2], that defines the fidelity
            level used when simulating recirculation.
                1: Recirc 1 pebble per recirculation loop, at a frequency
                defined by recirc_hz.  Simulation time and time step are
                manually updated to match 'real' time.
                2: Recirculate a group of pebbles each recirculation loop, 
                until recirc_target pebbles have been recirculated.
                This level also has an accelerated recirculation rate, and as
                such, the time steps in the simulation will not match
                'real' time.
        recirc_target : int
            For fidelity = 2 simulations, recirc_target is the minimum number
            of pebbles that must be recirculated.
        recirc_hz : float
            For fidelity = 1 simulations.  The frequency of pebble
            recirculation used, in units of [pebbles/second].
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
        self.fidelity = fidelity
        self.recirc_target = recirc_target
        self.recirc_hz = recirc_hz
        self.core_intake = core_intake
        self.core_main = core_main
        self.core_outtake = core_outtake
        self.recirc = recirc
        self.k_rate = k_rate
        self.down_flow = down_flow
        self.seed = seed
