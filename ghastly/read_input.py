import json
import numpy as np
from ghastly import core
from ghastly import simulation

rng = np.random.default_rng()


class InputBlock:
    '''
    Class for reading data from a JSON ghastly input file.
    '''

    def __init__(self, input_file):
        '''
        Initializes ghastly InputBlock object from data read from input_file.

        Parameters
        ----------
        input_file : JSON
            JSON file containing information necessary to run ghastly, see
            examples directory for example input files.  With the exception
            of core element names, such as "main_cyl" in the example input,
            variable names must match those in sample input.

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
        Using the InputBlock, create other ghastly objects needed to run
        LAMMPS and/or OpenMC.

        Returns
        -------
        sim_block: Sim object
            Ghastly Sim object with simulation-specific parameters and core
            zone dictionaries.

        '''

        core_intake = self.create_core_zone(self.core_intake_var)
        core_main = self.create_core_zone(self.core_main_var)
        core_outtake = self.create_core_zone(self.core_outtake_var)

        sim_block = self.create_sim_block(core_intake, core_main, core_outtake)

        return sim_block

    def create_core_zone(self, core_zone):
        '''
        Given a specific core_zone, such as intake, main, or outtake, create
        a dictionary with key:value pairs where each key is the name of a
        core element, and each value is the corresponding ghastly Core class
        object.

        Parameters
        ----------
        core_zone : dict
            Dictionary where each key:value pair corresponds to a single core
            element within the core_zone.  Values are also dictionaries with
            each key:value pair containing a parameter for the element and its
            value.

        Returns
        -------
        core_block : dict
            Dicitonary with key:value pairs where the key is the core element
            name, and the value is the corresponding ghastly Core object.

        '''

        core_block = {}
        for key, val in core_zone.items():
            if val["type"].casefold() == "cylinder":
                if val["open_bottom"] == True:
                    core_block[key] = core.CylCore(x_c=val["x_c"],
                                                   y_c=val["y_c"],
                                                   z_max=val["z_max"],
                                                   z_min=val["z_min"],
                                                   r=val["r"])
                else:
                    core_block[key] = core.CylCore(x_c=val["x_c"],
                                                   y_c=val["y_c"],
                                                   z_max=val["z_max"],
                                                   z_min=val["z_min"],
                                                   r=val["r"],
                                                   open_bottom="")

            elif val["type"].casefold() == "cone":
                if val["open_bottom"] == True:
                    core_block[key] = core.ConeCore(x_c=val["x_c"],
                                                    y_c=val["y_c"],
                                                    z_max=val["z_max"],
                                                    z_min=val["z_min"],
                                                    r_upper=val["r_upper"],
                                                    r_lower=val["r_lower"])
                else:
                    core_block[key] = core.ConeCore(x_c=val["x_c"],
                                                    y_c=val["y_c"],
                                                    z_max=val["z_max"],
                                                    z_min=val["z_min"],
                                                    r_upper=val["r_upper"],
                                                    r_lower=val["r_lower"],
                                                    open_bottom="")
            else:
                raise NameError("Type must be cylinder or cone.")

        return core_block

    def create_sim_block(self, core_intake, core_main, core_outtake):
        '''
        Creates a Sim object, using the intake, main, and outtake core zone
        dictionaries.  Default values for k_rate, down_flow, and seed are
        0.001, True, and a random integer between 1,000,000 and 100,000,000,
        respectively.

        Returns
        -------
        sim_block : Sim object
            Ghastly sim object containing simulation parameters and core
            objects.
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
            raise TypeError('''down_flow should be True for downward flow
                            and False for upward flow.''')

        match flow_case:
            case None:
                down_flow = True
            case _:
                down_flow = self.sim_var["down_flow"]

        seed_case = self.sim_var.get("seed")
        if type(seed_case) != int and seed_case != None:
            raise TypeError('''The random seed should be a large integer value,
                            or not in input block to randomly generate one.''')

        match seed_case:
            case None:
                seed = rng.integers(1000000, 100000000)
            case _:
                seed = self.sim_var["seed"]

        sim_block = simulation.Sim(r_pebble=self.sim_var["r_pebble"],
                                   t_final=self.sim_var["t_final"],
                                   pf=self.sim_var["pf"],
                                   core_intake=core_intake,
                                   core_main=core_main,
                                   core_outtake=core_outtake,
                                   k_rate=k_rate,
                                   down_flow=down_flow,
                                   seed=seed)
        return sim_block
