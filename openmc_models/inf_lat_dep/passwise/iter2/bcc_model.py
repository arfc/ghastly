import openmc
import openmc.deplete
import numpy as np


#materials
#graphite based on a3-3, triso layers pulled from reported values in
#Neutronics characteristics of a 165 MWth Xe-100 reactor, Mulder et al

#d_steps are (in days):
#d_steps = [1] + [4] + [4] + [10]*9 + [25]*10 + [50]*24
materials = openmc.Materials.from_xml('../materials.xml')
dep_file = 'i1-dep-res.h5'
res = openmc.deplete.Results(dep_file)
dep_t = res.get_times()

f = open('../nuc_list.txt', 'r')
nuc_w_data = []
for n in f.readlines():
    nuc_w_data.append(n.strip())

nuc_in_res = list(res[0].index_nuc.keys())

nuclides = list(set(nuc_w_data) & set(nuc_in_res))

print(nuclides)

p_bins = [(0,19),(19,26),(26,31),(31,36),(36,41),(41,None)]
p_comps = {'p01':{},
           'p12':{},
           'p23':{},
           'p34':{},
           'p45':{},
           'p56':{}}
p_t = [sum(dep_t[0:19]), 
       sum(dep_t[19:26]), 
       sum(dep_t[26:31]), 
       sum(dep_t[31:36]),
       sum(dep_t[36:41]), 
       sum(dep_t[41:])]
m_uco = materials[0].volume*materials[0].density

for nuc in nuclides:
    time, mass = res.get_mass('1', nuc, mass_units='g', time_units='d')
    for i, k in enumerate(list(p_comps.keys())):
        j = p_bins[i]
        m = sum(mass[j[0]:j[1]]*(time[j[0]:j[1]]/p_t[i]))/m_uco
        if m > 0:
            p_comps[k][nuc] = m
        else:
            pass

uco01= openmc.Material(name='UCO_01', material_id=10)
uco01.set_density('g/cm3', 10.4)
uco01.add_components(p_comps['p01'], percent_type = 'wo')
uco01.add_s_alpha_beta('c_Graphite')
uco01.depletable = False
uco01.temperature = 1159.15 #K

uco12= openmc.Material(name='UCO_12', material_id=20)
uco12.set_density('g/cm3', 10.4)
uco12.add_components(p_comps['p12'], percent_type = 'wo')
uco12.add_s_alpha_beta('c_Graphite')
uco12.depletable = False
uco12.temperature = 1159.15 #K

uco23= openmc.Material(name='UCO_23', material_id=30)
uco23.set_density('g/cm3', 10.4)
uco23.add_components(p_comps['p23'], percent_type = 'wo')
uco23.add_s_alpha_beta('c_Graphite')
uco23.depletable = False
uco23.temperature = 1159.15 #K

uco34= openmc.Material(name='UCO_34', material_id=40)
uco34.set_density('g/cm3', 10.4)
uco34.add_components(p_comps['p34'], percent_type = 'wo')
uco34.add_s_alpha_beta('c_Graphite')
uco34.depletable = False
uco34.temperature = 1159.15 #K

uco45= openmc.Material(name='UCO_45', material_id=50)
uco45.set_density('g/cm3', 10.4)
uco45.add_components(p_comps['p45'], percent_type = 'wo')
uco45.add_s_alpha_beta('c_Graphite')
uco45.depletable = False
uco45.temperature = 1159.15 #K

uco56= openmc.Material(name='UCO_56', material_id=60)
uco56.set_density('g/cm3', 10.4)
uco56.add_components(p_comps['p56'], percent_type = 'wo')
uco56.add_s_alpha_beta('c_Graphite')
uco56.depletable = False
uco56.temperature = 1159.15 #K

materials.append(uco01)
materials.append(uco12)
materials.append(uco23)
materials.append(uco34)
materials.append(uco45)
materials.append(uco56)

geometry = openmc.Geometry.from_xml('../geometry.xml', materials)

settings = openmc.Settings.from_xml('../settings.xml')

model = openmc.model.Model(geometry, materials, settings)
model.export_to_model_xml()

operator = openmc.deplete.CoupledOperator(model)

# 1549 effective full power days total lifetime
# dt steps based on table 1 in Adaptive Burnup..., Walter and Manera
# 165/1549 -> 1 EFPD ~= 0.1 MWd/kgHM
d_steps = [1] + [4] + [4] + [10]*9 + [25]*10 + [50]*24

reactor_power = 165.0*(10**6) #165 MWth, converted to W
#220K pebs, uco mass in 1 peb, wt percent of u in uco
uco_weight = 223000*m_uco*0.8945
specific_power = reactor_power/uco_weight #W/gHM

celi = openmc.deplete.CELIIntegrator(
        operator, d_steps, power_density=specific_power, timestep_units = 'd')

celi.integrate()


