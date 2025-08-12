import openmc
import openmc.deplete
import numpy as np

dep_file = 'i4-dep-res.h5'
mat_file = 'i4-mats.xml'
res = openmc.deplete.Results(dep_file)
dep_t = res.get_times()
step_comps = [res.export_to_materials(i, 
                                      path=mat_file)[0].get_nuclide_densities()
              for i in range(len(dep_t))]

#materials
#graphite based on a3-3, triso layers pulled from reported values in
#Neutronics characteristics of a 165 MWth Xe-100 reactor, Mulder et al

pass01 = {}
tot_time01 = sum(dep_t[0:19])
for i, step in enumerate(step_comps[0:19]):
    for k, v in step.items():
        if k in pass01:
            pass01[k] += v[1]*(dep_t[i]/tot_time01)
        else:
            pass01[k] = {}
            pass01[k] = v[1]*(dep_t[i]/tot_time01)

pass12 = {}
tot_time12 = sum(dep_t[19:26])
for i, step in enumerate(step_comps[0:19]):
    for k, v in step.items():
        if k in pass12:
            pass12[k] += v[1]*(dep_t[i+19]/tot_time12)
        else:
            pass12[k] = {}
            pass12[k] = v[1]*(dep_t[i+19]/tot_time12)

pass23 = {}
tot_time23 = sum(dep_t[26:31])
for i, step in enumerate(step_comps[26:31]):
    for k, v in step.items():
        if k in pass23:
            pass23[k] += v[1]*(dep_t[i+26]/tot_time23)
        else:
            pass23[k] = {}
            pass23[k] = v[1]*(dep_t[i+26]/tot_time23)

pass34 = {}
tot_time34 = sum(dep_t[31:36])
for i, step in enumerate(step_comps[31:36]):
    for k, v in step.items():
        if k in pass34:
            pass34[k] += v[1]*(dep_t[i+31]/tot_time34)
        else:
            pass34[k] = {}
            pass34[k] = v[1]*(dep_t[i+31]/tot_time34)

pass45 = {}
tot_time45 = sum(dep_t[36:41])
for i, step in enumerate(step_comps[36:41]):
    for k, v in step.items():
        if k in pass45:
            pass45[k] += v[1]*(dep_t[i+36]/tot_time45)
        else:
            pass45[k] = {}
            pass45[k] = v[1]*(dep_t[i+36]/tot_time45)

pass56 = {}
tot_time56 = sum(dep_t[41:])
for i, step in enumerate(step_comps[41:]):
    for k, v in step.items():
        if k in pass56:
            pass56[k] += v[1]*(dep_t[i+41]/tot_time56)
        else:
            pass56[k] = {}
            pass56[k] = v[1]*(dep_t[i+41]/tot_time56)


uco01= openmc.Material(name='UCO_01', material_id=7)
uco01.set_density('g/cm3', 10.4)
uco01.add_components(pass01, percent_type = 'ao')
uco01.add_s_alpha_beta('c_Graphite')
uco01.depletable = False
uco01.temperature = 1159.15 

uco12= openmc.Material(name='UCO_12', material_id=8)
uco12.set_density('g/cm3', 10.4)
uco12.add_components(pass12, percent_type = 'ao')
uco12.add_s_alpha_beta('c_Graphite')
uco12.depletable = False
uco12.temperature = 1159.15 

uco23= openmc.Material(name='UCO_23', material_id=9)
uco23.set_density('g/cm3', 10.4)
uco23.add_components(pass23, percent_type = 'ao')
uco23.add_s_alpha_beta('c_Graphite')
uco23.depletable = False
uco23.temperature = 1159.15 

uco34= openmc.Material(name='UCO_34', material_id=10)
uco34.set_density('g/cm3', 10.4)
uco34.add_components(pass34, percent_type = 'ao')
uco34.add_s_alpha_beta('c_Graphite')
uco34.depletable = False
uco34.temperature = 1159.15 

uco45= openmc.Material(name='UCO_45', material_id=11)
uco45.set_density('g/cm3', 10.4)
uco45.add_components(pass45, percent_type = 'ao')
uco45.add_s_alpha_beta('c_Graphite')
uco45.depletable = False
uco45.temperature = 1159.15 

uco56= openmc.Material(name='UCO_56', material_id=12)
uco56.set_density('g/cm3', 10.4)
uco56.add_components(pass56, percent_type = 'ao')
uco56.add_s_alpha_beta('c_Graphite')
uco56.depletable = False
uco56.temperature = 1159.15 

buffer = openmc.Material(name='buffer', material_id=14)
buffer.set_density('g/cm3', 1.05)
buffer.add_element('C', 1.0, percent_type='ao')
buffer.add_s_alpha_beta('c_Graphite')
buffer.temperature = 1159.15 #K

pyc = openmc.Material(name='PyC', material_id=15)
pyc.set_density('g/cm3', 1.9)
pyc.add_element('C', 1.0, percent_type='ao')
pyc.add_s_alpha_beta('c_Graphite')
pyc.temperature = 1159.15 #K

sic = openmc.Material(name='SiC', material_id=16)
sic.set_density('g/cm3', 3.2)
sic.add_element('C', 0.5, percent_type='ao')
sic.add_element('Si', 0.5, percent_type='ao')
sic.add_s_alpha_beta('c_Graphite')

graphite = openmc.Material(name='graphite', material_id=17)
graphite.set_density('kg/m3', 1700)
graphite.add_element('C', 1.0, percent_type='ao')
graphite.add_s_alpha_beta('c_Graphite')

he = openmc.Material(name='He', material_id=18)
he.set_density('atom/b-cm', 0.0006)
he.add_element('He', 1.0, percent_type='ao')
he.temperature = 778.15 #K

bcc_l = 100/(2697**(1/3)) #based on pf = 5394 pebs/m3

materials = openmc.Materials([uco01, uco12, uco23, uco34, uco45, uco56, 
                              buffer, pyc, sic, graphite, he])
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.verbosity = 7
settings.particles = 10000
settings.generations_per_batch = 10
settings.batches = 75
settings.inactive = 25
settings.seed = 463913357
settings.temperature = {'method' : 'interpolation', 'tolerance' : 10.0}
settings.export_to_xml()

tallies = openmc.Tallies()
mesh = openmc.RegularMesh()
mesh.dimension = [3, 3, 3]
mesh.lower_left = [-0.5*bcc_l, -0.5*bcc_l, -0.5*bcc_l]
mesh.upper_right = [0.5*bcc_l, 0.5*bcc_l, 0.5*bcc_l]
mesh_filter = openmc.MeshFilter(mesh)
spatial = openmc.Tally(name='spatial')
spatial.filters = [mesh_filter]
spatial.scores = ['flux', 'fission']
heating = openmc.Tally(name='heating')
heating.scores = ['heating']
spectrum = openmc.Tally(name='spectrum')
ene = np.logspace(np.log10(1e-5), np.log10(20.0e6), 501)
ene_filter = openmc.EnergyFilter(ene)
spectrum.filters = [ene_filter]
spectrum.scores = ['flux']
tallies.append(spatial)
tallies.append(heating)
tallies.append(spectrum)
tallies.export_to_xml()
