class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''

    def __init__(self, r_pebble, t_final, down_flow=True):
        '''
        Initializes the Sim class.
        '''
        self.r_pebble = r_pebble
        self.t_final = t_final
        self.down_flow = down_flow
