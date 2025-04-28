import numpy as np
import openmc
import openmc.deplete


class Mat():
    '''
    Class containing OpenMC material metadata necessary for pebble tracking.
    '''

    def __init__(self, material, pass_num, reg_id, l_type=0):
        '''
        Initializes a single instance of the Mat class.

        Parameters
        ----------
        material : OpenMC Material
            The OpenMC material object this metadata applies to.
        pass_num : int
            Integer corresponding to the number of passes a pebble with this
            material would have completed.  A fresh pebble would have a
            pass_num of 0, for example.
        reg_id  : str
            Region identifier, `reg_id`, for the region this material 
            corresponds to.
        l_type : int
            The type in LAMMPS that corresponds to this material.  These will
            be determined by ghastly during simulation set up.
        '''
        self.material = material
        self.pass_num = pass_num
        self.reg_id = reg_id
        self.l_type = l_type

def avg_comp(dep_file, mat_file):
    '''
    given a depletion results h5 file and mat file path, returns
    a dict containing the burnup-step-wise average compositions for
    depleted fuel material, in the same units as the dep results file.
    currently assumes first mat is the fuel - look for fix so function
    can find which one is the fuel on its own later
    '''

    res = openmc.deplete.Results(dep_file)
    n = len(res.get_times())
    comps = [res.export_to_materials(i,
                                     path=mat_file)[0].get_nuclide_densities() 
             for i in range(n)]
    avg_comp = {}
    for comp in comps:
        for k, v in comp.items():
            if k in avg_comp:
                avg_comp[k]['iso'] += v[1]
                avg_comp[k]['count'] += 1

            else:
                avg_comp[k] = {}
                avg_comp[k]['iso'] = v[1]
                avg_comp[k]['count'] = 1
    core_avg = {}
    for k, v in avg_comp.items():
        core_avg[k] = v['iso']/v['count']

    return core_avg
