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

pass_ids = [01, 12, 23, 34, 45, 56]
dep_i = [(0, 19),
         (19, 26),
         (26, 31),
         (31, 36),
         (36, 41)]
passes = {}
tot_times = [sum(dep_t[i[0]:i[1]]) for i in dep_i]

for i, interval in enumerate(dep_i):
    comp = {}
    for j, step in enumerate(step_comps[interval[0]:interval[1]]):
        for k, v in step.items():
            if k in comp:
                comp[k] += v[1]*(dep_t[j]/tot_time[i])
            else:
                comp[k] = v[1]*(dep_t[j]/tot_time[i])
    passes[pass_ids[i]] = comp
ucos = []
uid = 7
for p_id in pass_ids:
    mat_name = 'UCO_'+str(p_id)
    uco= openmc.Material(name=mat_name, material_id = uid)
    uco.set_density('g/cm3', 10.4)
    uco.add_components(passes[p_id], percent_type = 'ao')
    uco.add_s_alpha_beta('c_Graphite')
    uco.depletable = False
    uco.temperature = 1159.15
    ucos.append(uco)
    uid += 1

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

materials = openmc.Materials([ucos, buffer, pyc, sic, graphite, he])
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.verbosity = 7
settings.particles = 20000
settings.batches = 80
settings.inactive = 40
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
