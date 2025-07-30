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
#also doublecheck isotopic compostions conventions in mats
#graphite based on a3-3, triso layers pulled from reported values in
#Neutronics characteristics of a 165 MWth Xe-100 reactor, Mulder et al

#from below, d_steps are as follows (in [days]), so the pass-wise
#indices would be (remember open on end, each pass ~= 258 days
#pass1 : 0, then everything from 1 to 10 days, inclusive, then the first 6
#steps of [25] days : 
#d_steps = [1] + [4] + [4] + [10]*9 + [25]*10 + [50]*24


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

uco01= openmc.Material(name='UCO_01')
uco01.set_density('g/cm3', 10.4)
uco01.add_components(pass01, percent_type = 'ao')

uco12= openmc.Material(name='UCO_12')
uco12.set_density('g/cm3', 10.4)
uco12.add_components(pass12, percent_type = 'ao')

uco23= openmc.Material(name='UCO_23')
uco23.set_density('g/cm3', 10.4)
uco23.add_components(pass23, percent_type = 'ao')

uco34= openmc.Material(name='UCO_34')
uco34.set_density('g/cm3', 10.4)
uco34.add_components(pass34, percent_type = 'ao')

uco45= openmc.Material(name='UCO_45')
uco45.set_density('g/cm3', 10.4)
uco45.add_components(pass45, percent_type = 'ao')

uco56= openmc.Material(name='UCO_56')
uco56.set_density('g/cm3', 10.4)
uco56.add_components(pass56, percent_type = 'ao')

ucoavg = openmc.Material.mix_materials([uco01, uco12,
                                        uco23, uco34,
                                        uco45, uco56], [0.16, 0.17,
                                                          0.17, 0.17,
                                                          0.17, 0.16], 'vo',
                                       material_id=13)
ucoavg.add_s_alpha_beta('c_Graphite')
ucoavg.depletable = False
ucoavg.temperature = 1159.15 #K

comp1549days = {}
for k, v in step_comps[-1].items():
        comp1549days[k] = {}
        comp1549days[k] = v[1]

ucot = openmc.Material(name='UCO_TRACKED', material_id=14)
ucot.set_density('g/cm3', 10.4)
ucot.add_components(comp1549days, percent_type = 'ao')
ucot.add_s_alpha_beta('c_Graphite')
ucot.depletable = True
ucot.temperature = 1159.15 #K

buffer = openmc.Material(name='buffer', material_id=15)
buffer.set_density('g/cm3', 1.05)
buffer.add_element('C', 1.0, percent_type='ao')
buffer.add_s_alpha_beta('c_Graphite')
buffer.temperature = 1159.15 #K

pyc = openmc.Material(name='PyC', material_id=16)
pyc.set_density('g/cm3', 1.9)
pyc.add_element('C', 1.0, percent_type='ao')
pyc.add_s_alpha_beta('c_Graphite')
pyc.temperature = 1159.15 #K

sic = openmc.Material(name='SiC', material_id=17)
sic.set_density('g/cm3', 3.2)
sic.add_element('C', 0.5, percent_type='ao')
sic.add_element('Si', 0.5, percent_type='ao')
sic.add_s_alpha_beta('c_Graphite')

graphite = openmc.Material(name='graphite', material_id=18)
graphite.set_density('kg/m3', 1700)
graphite.add_element('C', 1.0, percent_type='ao')
graphite.add_s_alpha_beta('c_Graphite')

he = openmc.Material(name='He', material_id=19)
he.set_density('atom/b-cm', 0.0006)
he.add_element('He', 1.0, percent_type='ao')
he.temperature = 778.15 #K



#geometry parameters - remember LAMMPS coords are in meters, this is in cm
#remember, triso order: kernel, buffer, pyc, sic, pyc
peb_or = 3.0 #outer radius of whole pebble
peb_ir = 2.5 # radius of the region that has trisos in it only
triso_r = [0.02125, 0.03125, 0.03525, 0.03875, 0.04275]
bcc_l = 100/(2697**(1/3)) #based on pf = 5394 pebs/m3

materials = openmc.Materials([ucoavg, ucot, buffer, pyc, sic, graphite, he])
openmc.Materials(materials).export_to_xml()

geometry = openmc.Geometry.from_xml("geometry.xml", materials)

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.verbosity = 7
settings.particles = 5000
settings.generations_per_batch = 7
settings.batches = 100
settings.inactive = 70
settings.seed = 463913357
settings.temperature = {'method' : 'interpolation', 'tolerance' : 10.0}
#settings.volume_calculations = [vol_calc]
settings.export_to_xml()

tallies = openmc.Tallies()
mesh = openmc.RegularMesh()
mesh.dimension = [2,2,2]
mesh.lower_left = [-0.5*bcc_l, -0.5*bcc_l, -0.5*bcc_l]
mesh.upper_right = [0.5*bcc_l, 0.5*bcc_l, 0.5*bcc_l]
mesh_filter = openmc.MeshFilter(mesh)
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux', 'fission']
tallies.append(tally)
tallies.export_to_xml()

openmc.run()
