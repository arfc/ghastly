import openmc
import openmc.deplete
import numpy as np

dep_file = '../i4-dep-res.h5'
mat_file = '../i4-mats.xml'
res = openmc.deplete.Results(dep_file)
dep_t = res.get_times()
step_comps = [res.export_to_materials(i, 
                                      path=mat_file)[0].get_nuclide_densities()
              for i in range(len(dep_t))]

comp999days = {}
for k, v in step_comps[35].items():
        comp999days[k] = {}
        comp999days[k] = v[1]

ucot = openmc.Material(name='UCO_TRACKED', material_id=13)
ucot.set_density('g/cm3', 10.4)
ucot.add_components(comp999days, percent_type = 'ao')
ucot.add_s_alpha_beta('c_Graphite')
ucot.depletable = True
ucot.temperature = 1159.15 #K

materials = openmc.Materials.from_xml('../materials.xml')
materials.append(ucot)
geometry = openmc.Geometry.from_xml("../geometry.xml", materials)
settings = openmc.Settings.from_xml('../settings.xml')
tallies = openmc.Tallies.from_xml('../tallies.xml')

model = openmc.model.Model(geometry, materials, settings, tallies)
model.export_to_model_xml()

openmc.run()
