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
        if "recirc" in params:
            self.recirc_var = params["recirc"]
        else:
            self.recirc_var = {}
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
                core_block[key] = core.CylCore(x_c=val["x_c"],
                                               y_c=val["y_c"],
                                               z_max=val["z_max"],
                                               z_min=val["z_min"],
                                               r=val["r"],
                                               open_bottom=val["open_bottom"],
                                               open_top=val["open_top"])

            elif val["type"].casefold() == "cone":
                core_block[key] = core.ConeCore(x_c=val["x_c"],
                                                y_c=val["y_c"],
                                                z_max=val["z_max"],
                                                z_min=val["z_min"],
                                                r_upper=val["r_upper"],
                                                r_lower=val["r_lower"],
                                                open_bottom=val["open_bottom"],
                                                open_top=val["open_top"])
            else:
                raise NameError("Type must be cylinder or cone.")

        return core_block

    def create_recirc_zone(self, recirc_zone):
        '''
        Given the recirc block from a ghastly input file, create
        a dictionary with key:value pairs where each key is the name of a
        recirc element, and each value is the corresponding ghastly Core class
        object.

        Parameters
        ----------
        recirc_zone : dict
            Dictionary where each key:value pair corresponds to a single
            recirc element within the core_zone, generally corresponding to an
            inlet and outlet.  Values are also dictionaries with
            each key:value pair containing the name of a region parameter 
            and its value.

        Returns
        -------
        recirc_block : dict
            Dicitonary with key:value pairs where the key is the recirc element
            name, and the value is the corresponding ghastly Core object.

        '''

        recirc_block = {}
        for key, val in recirc_zone.items():
            if val["type"].casefold() == "cylinder":
                recirc_block[key] = core.CylCore(x_c=val["x_c"],
                                               y_c=val["y_c"],
                                               z_max=val["z_max"],
                                               z_min=val["z_min"],
                                               r=val["r"])

            elif val["type"].casefold() == "cone":
                recirc_block[key] = core.ConeCore(x_c=val["x_c"],
                                               y_c=val["y_c"],
                                               z_max=val["z_max"],
                                               z_min=val["z_min"],
                                               r_upper=val["r_upper"],
                                                 r_lower=val["r_lower"])
            else:
                raise NameError("Type must be cylinder or cone.")

        return recirc_block


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

        hz_case = self.sim_var.get("recirc_hz")
        match hz_case:
            case None:
                recirc_hz = 1
            case _:
                recirc_hz = self.sim_var["recirc_hz"]
        
        target_case = self.sim_var.get("recirc_target")
        if type(target_case) != int and target_case != None:
            raise TypeError('''The target number of pebbles to recirculate
                            should be an integer''')
        match target_case:
            case None:
                recirc_target = 1
            case _:
                recirc_target = self.sim_var["recirc_target"]

        fidelity_case = self.sim_var.get("fidelity")
        if type(fidelity_case) != int and fidelity_case != None:
            raise TypeError('''The fidelity level must be one of the following
                            integers: [1, 2]''')
        match fidelity_case:
            case None:
                fidelity = 1
            case _:
                fidelity = self.sim_var["fidelity"]


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
                seed = rng.integers(1000000, 9999999)
            case _:
                seed = self.sim_var["seed"]

        match len(self.recirc_var):
            case 0:
                recirc = {}
            case _:
                recirc = self.create_recirc_zone(self.recirc_var)

        sim_block = simulation.Sim(r_pebble=self.sim_var["r_pebble"],
                                   t_final=self.sim_var["t_final"],
                                   pf=self.sim_var["pf"],
                                   fidelity=fidelity,
                                   recirc_target=recirc_target,
                                   recirc_hz=recirc_hz,
                                   core_intake=core_intake,
                                   core_main=core_main,
                                   core_outtake=core_outtake,
                                   recirc=recirc,
                                   k_rate=k_rate,
                                   down_flow=down_flow,
                                   seed=seed)
        return sim_block
