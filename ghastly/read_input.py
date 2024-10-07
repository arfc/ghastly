import json
from ghastly import core
from ghastly import simulation

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

    def create_core_zone(self, core_zone):
        '''
        sub-function in create_obj to make the core objects specfically
        '''

        core_block = {}
        for key, val in core_zone.items():
            if val["type"] == "cylinder":
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
                    
            elif val["type"] == "cone":
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

    def create_obj(self):
        '''
        using input file, create ghastly objects need to run simulation
        '''
        pass
        #start with core parts, get core main/outtake/intake dicts, so then
        #you can pass them to sim class.

        core_intake = create_core_zone(self.core_intake_var)
        core_main = create_core_zone(self.core_main_var)
        core_outtake = create_core_zone(self.core_outtake_var)

        sim_params = simulation.Sim(r_pebble = self.sim_var["r_pebble"],
                                    t_final = self.sim_var["t_final"],
                                    pf = self.sim_var["pf"],
                                    core_intake = core_intake,
                                    core_main = core_main,
                                    core_outtake = core_outtake,
                                    k_rate = self.sim_var.get("k_rate"),
                                    down_flow = self.sim_var.get("down_flow"),
                                    seed = self.seed.get("seed"))

        return sim_params


        
