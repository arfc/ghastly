class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''

    def __init__(self, r_pebble, t_final, down_flow=True):
        '''
        Initializes the Sim class.

        Parameters
        ----------
        r_pebble : float
            Radius of the pebbles in the simulation [m].
        t_final : float
            Total reactor-time that is being simulated [s].
        down_flow : bool
            The direction of flow in the reactor core.  True means it is
            flowing downward, False means it is flowing upward.  Default is
            down.
        '''
        self.r_pebble = r_pebble
        self.t_final = t_final
        self.down_flow = down_flow
