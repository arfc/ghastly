import json
import numpy as np
from ghastly import core
from ghastly import simulation

rng = np.random.default_rng()

class InputBlock:
    '''
    class for reading data from a ghastly input file (json)
    '''
    def __init__(self, input_file):
        '''
        initializes ghastly input object using a ghastly input file.
        '''

        f = open(input_file, 'r')
        params = json.load(f)
        
        self.sim_var = params["simulation"]
        self.core_intake_var = params["core_intake"]
        self.core_main_var = params["core_main"]
        self.core_outtake_var = params["core_outtake"]
        self.lammps_var = params["lammps_var"]


    def create_obj(self):
        '''
        using input file, create ghastly objects need to run simulation
        '''
        #start with core parts, get core main/outtake/intake dicts, so then
        #you can pass them to sim class.

        core_intake = self.create_core_zone(self.core_intake_var)
        core_main = self.create_core_zone(self.core_main_var)
        core_outtake = self.create_core_zone(self.core_outtake_var)

        sim_block = self.create_sim_block(core_intake, core_main, core_outtake)
        
        return sim_block

    def create_core_zone(self, core_zone):
        '''
        sub-function in create_obj to make the core objects specfically
        '''

        core_block = {}
        for key, val in core_zone.items():
            if val["type"].casefold() == "cylinder":
                if val["open_bottom"] == True:
                    core_block[key] = core.CylCore(x_c = val["x_c"],
                                                   y_c = val["y_c"],
                                                   z_max = val["z_max"],
                                                   z_min = val["z_min"],
                                                   r = val["r"])
                else:
                    core_block[key] = core.CylCore(x_c = val["x_c"],
                                                   y_c = val["y_c"],
                                                   z_max = val["z_max"],
                                                   z_min = val["z_min"],
                                                   r = val["r"],
                                                   open_bottom = "")
                    
            elif val["type"].casefold() == "cone":
                if val["open_bottom"] == True:
                    core_block[key] = core.ConeCore(x_c = val["x_c"],
                                                    y_c = val["y_c"],
                                                    z_max = val["z_max"],
                                                    z_min = val["z_min"],
                                                    r_upper = val["r_upper"],
                                                    r_lower = val["r_lower"])
                else:
                    core_block[key] = core.ConeCore(x_c = val["x_c"],
                                                    y_c = val["y_c"],
                                                    z_max = val["z_max"],
                                                    z_min = val["z_min"],
                                                    r_upper = val["r_upper"],
                                                    r_lower = val["r_lower"],
                                                     open_bottom = "")
            else:
                raise NameError("Type must be cylinder or cone.")

        return core_block
        
    def create_sim_block(self, core_intake, core_main, core_outtake):
        '''
        create sim class object
        '''
        k_case = self.sim_var.get("k_rate")
        if type(k_case) != float and k_case != None:
            raise TypeError('''The contraction rate should be a value between 0
                            and 1, non-inclusive.''')
        match k_case:
            case None:
                k_rate = 0.001
            case _:
                k_rate = self.sim_var["k_rate"]

        flow_case = self.sim_var.get("down_flow")
        if type(flow_case) != bool and flow_case != None:
            raise TypeError('''down_flow should be true for downward flow
                            and false for upward flow''')

        match flow_case:
            case None:
                down_flow = True
            case _:
                down_flow = self.sim_var["down_flow"]

        seed_case = self.sim_var.get("seed")
        if type(seed_case) != int and seed_case != None:
            raise TypeError('''The random seed should be a large integer value,
                            or not in input block to randomly generate one''')

        match seed_case:
            case None:
                seed = rng.integers(1000000,100000000)
            case _:
                seed = self.sim_var["seed"]

        sim_block = simulation.Sim(r_pebble = self.sim_var["r_pebble"],
                                   t_final = self.sim_var["t_final"],
                                   pf = self.sim_var["pf"],
                                   core_intake = core_intake,
                                   core_main = core_main,
                                   core_outtake = core_outtake,
                                   k_rate = k_rate,
                                   down_flow = down_flow,
                                   seed = seed)
        return sim_block


