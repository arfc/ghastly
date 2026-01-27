import openmc
import openmc.deplete
import numpy as np

ucot = openmc.Material(name='UCO_TRACKED', material_id=14)
ucot.set_density('g/cm3', 10.4)
ucot.add_nuclide("U235", 0.1386, percent_type='wo')
ucot.add_nuclide("U238",0.7559, percent_type='wo')
ucot.add_element("O", 0.06025, percent_type='wo')
ucot.add_element('C', 0.04523, percent_type='wo')
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
