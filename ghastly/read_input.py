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

    def create_obj(self):
        '''
        using input file, create ghastly objects need to run simulation
        '''

        #start with core parts, get core main/outtake/intake dicts, so then
        #you can pass them to sim class. 
        
        core_intake = {}
        core_main = {}
        core_outtake = {}
        for key, val in self.core_intake.items():
            if val["type"] == "cylinder":
                core_intake[key] = core.CylCore(x_c = val["x_c"],
                                                y_c = val["y_c"],
                                                z_max = val["z_max"],
                                                z_min = val["z_min"],
                                                r = val["r"])
            elif val["type"] == "cone":
                #ghastly conecore
                pass
            else:
                raise NameError("Type must be cylinder or cone.")
