import openmc
import openmc.deplete
import numpy as np

dep_file = 'i0-dep-res.h5'
mat_file = 'input_mats.xml'
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


#realized something - discuss most logical option with Luke??
# if I have "fresh" as a pass, and weight it like the other passes, it will
#skew the results towards fresh, bc fresh doesn't exist the whole time (just
#like how the last step of pass 5 -> 6, the most burnt, doesn't always exist, 
# "fresh" is just the first step of pass 0->1
temp_comp = {}
for step in step_comps[0:19]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass01 = {}
for k, v in temp_comp.items():
    pass01[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[19:26]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1      
pass12 = {}
for k, v in temp_comp.items():
    pass12[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[26:31]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass23 = {}
for k, v in temp_comp.items():
    pass23[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[31:36]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass34 = {}
for k, v in temp_comp.items():
    pass34[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[36:41]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass45 = {}
for k, v in temp_comp.items():
    pass45[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[41:]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass56 = {}
for k, v in temp_comp.items():
    pass56[k] = v['iso']/v['count']

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
                                                          0.17, 0.16], 'vo')
ucoavg.add_s_alpha_beta('c_Graphite')
ucoavg.depletable = False
ucoavg.temperature = 1159.15 #K



ucot = openmc.Material(name='UCO_TRACKED')
ucot.set_density('g/cm3', 10.4)
ucot.add_nuclide("U235", 0.1386, percent_type='wo')
ucot.add_nuclide("U238",0.7559, percent_type='wo')
ucot.add_element("O", 0.06025, percent_type='wo')
ucot.add_element('C', 0.04523, percent_type='wo')
ucot.add_s_alpha_beta('c_Graphite')
ucot.depletable = True
ucot.temperature = 1159.15 #K

buffer = openmc.Material(name='buffer')
buffer.set_density('g/cm3', 1.05)
buffer.add_element('C', 1.0, percent_type='ao')
buffer.add_s_alpha_beta('c_Graphite')
buffer.temperature = 1159.15 #K

pyc = openmc.Material(name='PyC')
pyc.set_density('g/cm3', 1.9)
pyc.add_element('C', 1.0, percent_type='ao')
pyc.add_s_alpha_beta('c_Graphite')
pyc.temperature = 1159.15 #K

sic = openmc.Material(name='SiC')
sic.set_density('g/cm3', 3.2)
sic.add_element('C', 0.5, percent_type='ao')
sic.add_element('Si', 0.5, percent_type='ao')
sic.add_s_alpha_beta('c_Graphite')

graphite = openmc.Material(name='graphite')
graphite.set_density('kg/m3', 1700)
graphite.add_element('C', 1.0, percent_type='ao')
graphite.add_s_alpha_beta('c_Graphite')

he = openmc.Material(name='He')
he.set_density('atom/b-cm', 0.0006)
he.add_element('He', 1.0, percent_type='ao')
he.temperature = 778.15 #K



#geometry parameters - remember LAMMPS coords are in meters, this is in cm
#remember, triso order: kernel, buffer, pyc, sic, pyc
peb_or = 3.0 #outer radius of whole pebble
peb_ir = 2.5 # radius of the region that has trisos in it only
triso_r = [0.02125, 0.03125, 0.03525, 0.03875, 0.04275]
bcc_l = 100/(2697**(1/3)) #based on pf = 5394 pebs/m3
c_coord = 0.5*bcc_l



#define the edges of the unit cell
bcc_min_x = openmc.XPlane(x0=-c_coord, boundary_type='reflective')
bcc_max_x = openmc.XPlane(x0=c_coord, boundary_type='reflective')
bcc_min_y = openmc.YPlane(y0=-c_coord, boundary_type='reflective')
bcc_max_y = openmc.YPlane(y0=c_coord, boundary_type='reflective')
bcc_min_z = openmc.ZPlane(z0=-c_coord, boundary_type='reflective')
bcc_max_z = openmc.ZPlane(z0=c_coord, boundary_type='reflective')

#define the cells and universe for the triso particles
triso_reg = [openmc.Sphere(r=r) for r in triso_r[:-1]]
tr_triso_cells = [openmc.Cell(fill=ucot, region=-triso_reg[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_reg[0] & -triso_reg[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_reg[1] & -triso_reg[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_reg[2] & -triso_reg[3]),
               openmc.Cell(fill=pyc, region=+triso_reg[3])]
tr_triso_univ = openmc.Universe(cells=tr_triso_cells)

bg_triso_cells = [openmc.Cell(fill=ucoavg, region=-triso_reg[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_reg[0] & -triso_reg[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_reg[1] & -triso_reg[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_reg[2] & -triso_reg[3]),
               openmc.Cell(fill=pyc, region=+triso_reg[3])]
bg_triso_univ = openmc.Universe(cells=bg_triso_cells)



# define where the particles will be packed
# BCC lattice has one sphere at the center, and a 1/4 sphere at each corner, 
#tight enough they are all touching. that part should be handled with your 
#calculation of the cube length (the body diameter = 2 peb diameters)

#start with the easy part: the center pebble

body_peb_in = openmc.Sphere(r = peb_ir)
body_wfuel_bound = -body_peb_in
body_peb_out = openmc.Sphere(r=peb_or)
body_nofuel_reg = -body_peb_out & +body_peb_in

#remember, this is the fueled part only, so it uses peb_ir.  
#Center the BCC at 0.
#generate the triso centers in the body pebble:

body_centers = openmc.model.pack_spheres(triso_r[-1], region=body_wfuel_bound,
                                         num_spheres=19000, seed = 978397880)
body_trisos = [openmc.model.TRISO(triso_r[-1], tr_triso_univ, center) 
               for center in body_centers]

#now define the no-fueled outer shell for the body peb

body_wfuel = openmc.Cell(region=body_wfuel_bound)
lower_left, upper_right = body_wfuel.region.bounding_box
shape = (10, 10, 10)
pitch = (upper_right - lower_left)/shape
body_lattice = openmc.model.create_triso_lattice(
    body_trisos, lower_left, pitch, shape, graphite)

body_wfuel.fill = body_lattice
body_nofuel = openmc.Cell(fill=graphite, region=body_nofuel_reg)

body_cells = [body_wfuel, body_nofuel]



#this will follow the same basic steps as the center pebble.  initial triso 
#definition shouldn't change.
#this one is the -1,-1,-1 corner
c1_peb_in = openmc.Sphere(x0 = -c_coord, 
                          y0 =-c_coord, 
                          z0 =-c_coord, 
                          r = peb_ir)
c1_peb_out = openmc.Sphere(x0 = -c_coord, 
                           y0 =-c_coord, 
                           z0 =-c_coord,
                           r = peb_or)
c1_corner_bound = +bcc_min_x & +bcc_min_y & +bcc_min_z
c1_wfuel_reg = -c1_peb_in & c1_corner_bound
c1_nofuel_reg = -c1_peb_out & +c1_peb_in & c1_corner_bound

c1_centers = openmc.model.pack_spheres(triso_r[-1], region=-c1_peb_in, 
                                       num_spheres=19000, seed = 987432789)
c1_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c1_centers]

c1_wfuel = openmc.Cell(region=c1_wfuel_reg)
lower_left, upper_right = c1_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c1_lattice = openmc.model.create_triso_lattice(
    c1_trisos, lower_left, pitch, shape, graphite)

c1_wfuel.fill = c1_lattice
c1_nofuel = openmc.Cell(fill=graphite, region=c1_nofuel_reg)

c1_cells = [c1_wfuel, c1_nofuel]



#c2: -1,-1,1
c2_peb_in = openmc.Sphere(x0 = -c_coord, y0 =-c_coord, z0 =c_coord, r = peb_ir)
c2_peb_out = openmc.Sphere(x0 = -c_coord, y0 =-c_coord, z0 =c_coord,r = peb_or)
c2_corner_bound = +bcc_min_x & +bcc_min_y & -bcc_max_z
c2_wfuel_reg = -c2_peb_in & c2_corner_bound
c2_nofuel_reg = -c2_peb_out & +c2_peb_in & c2_corner_bound

c2_centers = openmc.model.pack_spheres(triso_r[-1], region=-c2_peb_in, 
                                       num_spheres=19000, seed = 895490889)
c2_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c2_centers]

c2_wfuel = openmc.Cell(region=c2_wfuel_reg)
lower_left, upper_right = c2_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c2_lattice = openmc.model.create_triso_lattice(
    c2_trisos, lower_left, pitch, shape, graphite)

c2_wfuel.fill = c2_lattice
c2_nofuel = openmc.Cell(fill=graphite, region=c2_nofuel_reg)

c2_cells = [c2_wfuel, c2_nofuel]



#c3: -1, 1, -1
c3_peb_in = openmc.Sphere(x0 = -c_coord, 
                          y0 =c_coord, 
                          z0 =-c_coord,
                          r = peb_ir)
c3_peb_out = openmc.Sphere(x0 = -c_coord, 
                           y0 =c_coord, 
                           z0 =-c_coord,
                           r = peb_or)
c3_corner_bound = +bcc_min_x & -bcc_max_y & +bcc_min_z
c3_wfuel_reg = -c3_peb_in & c3_corner_bound
c3_nofuel_reg = -c3_peb_out & +c3_peb_in & c3_corner_bound

c3_centers = openmc.model.pack_spheres(triso_r[-1], region=-c3_peb_in, 
                                       num_spheres=19000, seed = 356657913)
c3_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c3_centers]

c3_wfuel = openmc.Cell(region=c3_wfuel_reg)
lower_left, upper_right = c3_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c3_lattice = openmc.model.create_triso_lattice(
    c3_trisos, lower_left, pitch, shape, graphite)

c3_wfuel.fill = c3_lattice
c3_nofuel = openmc.Cell(fill=graphite, region=c3_nofuel_reg)

c3_cells = [c3_wfuel, c3_nofuel]



#c4: -1, 1, 1
c4_peb_in = openmc.Sphere(x0 = -c_coord, 
                          y0 =c_coord,
                          z0 =c_coord, 
                          r = peb_ir)
c4_peb_out = openmc.Sphere(x0 = -c_coord, 
                           y0 =c_coord, 
                           z0 =c_coord,
                           r = peb_or)
c4_corner_bound = +bcc_min_x & -bcc_max_y & -bcc_max_z
c4_wfuel_reg = -c4_peb_in & c4_corner_bound
c4_nofuel_reg = -c4_peb_out & +c4_peb_in & c4_corner_bound


c4_centers = openmc.model.pack_spheres(triso_r[-1], region=-c4_peb_in, 
                                       num_spheres=19000, seed = 353623684)
c4_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c4_centers]

c4_wfuel = openmc.Cell(region=c4_wfuel_reg)
lower_left, upper_right = c4_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c4_lattice = openmc.model.create_triso_lattice(
    c4_trisos, lower_left, pitch, shape, graphite)

c4_wfuel.fill = c4_lattice
c4_nofuel = openmc.Cell(fill=graphite, region=c4_nofuel_reg)

c4_cells = [c4_wfuel, c4_nofuel]



#c5: 1, -1, -1
c5_peb_in = openmc.Sphere(x0 = c_coord, 
                          y0 =-c_coord, 
                          z0 =-c_coord, 
                          r = peb_ir)
c5_peb_out = openmc.Sphere(x0 = c_coord, 
                           y0 =-c_coord, 
                           z0 =-c_coord,
                           r = peb_or)
c5_corner_bound = -bcc_max_x & +bcc_min_y & +bcc_min_z
c5_wfuel_reg = -c5_peb_in & c5_corner_bound
c5_nofuel_reg = -c5_peb_out & +c5_peb_in & c5_corner_bound

c5_centers = openmc.model.pack_spheres(triso_r[-1], region=-c5_peb_in, 
                                       num_spheres=19000, seed = 897578459)
c5_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c5_centers]

c5_wfuel = openmc.Cell(region=c5_wfuel_reg)
lower_left, upper_right = c5_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c5_lattice = openmc.model.create_triso_lattice(
    c5_trisos, lower_left, pitch, shape, graphite)

c5_wfuel.fill = c5_lattice
c5_nofuel = openmc.Cell(fill=graphite, region=c5_nofuel_reg)

c5_cells = [c5_wfuel, c5_nofuel]



#c6: 1, -1, 1
c6_peb_in = openmc.Sphere(x0 = c_coord, 
                          y0 =-c_coord, 
                          z0 =c_coord, 
                          r = peb_ir)
c6_peb_out = openmc.Sphere(x0 = c_coord, 
                           y0 =-c_coord,
                           z0 =c_coord,
                           r = peb_or)
c6_corner_bound = -bcc_max_x & +bcc_min_y & -bcc_max_z
c6_wfuel_reg = -c6_peb_in & c6_corner_bound
c6_nofuel_reg = -c6_peb_out & +c6_peb_in & c6_corner_bound

c6_centers = openmc.model.pack_spheres(triso_r[-1], region=-c6_peb_in, 
                                       num_spheres=19000, seed = 429875621)
c6_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c6_centers]

c6_wfuel = openmc.Cell(region=c6_wfuel_reg)
lower_left, upper_right = c6_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c6_lattice = openmc.model.create_triso_lattice(
    c6_trisos, lower_left, pitch, shape, graphite)

c6_wfuel.fill = c6_lattice
c6_nofuel = openmc.Cell(fill=graphite, region=c6_nofuel_reg)

c6_cells = [c6_wfuel, c6_nofuel]

#c7: 1, 1, -1
c7_peb_in = openmc.Sphere(x0 = c_coord, 
                          y0 =c_coord, 
                          z0 =-c_coord, 
                          r = peb_ir)
c7_peb_out = openmc.Sphere(x0 = c_coord, 
                           y0 =c_coord, 
                           z0 =-c_coord,
                           r = peb_or)
c7_corner_bound = -bcc_max_x & -bcc_max_y & +bcc_min_z
c7_wfuel_reg = -c7_peb_in & c7_corner_bound
c7_nofuel_reg = -c7_peb_out & +c7_peb_in & c7_corner_bound

c7_centers = openmc.model.pack_spheres(triso_r[-1], region=-c7_peb_in, 
                                       num_spheres=19000, seed = 489638795)
c7_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c7_centers]

c7_wfuel = openmc.Cell(region=c7_wfuel_reg)
lower_left, upper_right = c7_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c7_lattice = openmc.model.create_triso_lattice(
    c7_trisos, lower_left, pitch, shape, graphite)

c7_wfuel.fill = c7_lattice
c7_nofuel = openmc.Cell(fill=graphite, region=c7_nofuel_reg)

c7_cells = [c7_wfuel, c7_nofuel]

#c8: 1, 1, 1
c8_peb_in = openmc.Sphere(x0 = c_coord, 
                          y0 =c_coord, 
                          z0 =c_coord, 
                          r = peb_ir)
c8_peb_out = openmc.Sphere(x0 = c_coord, 
                           y0 =c_coord, 
                           z0 =c_coord,
                           r = peb_or)
c8_corner_bound = -bcc_max_x & -bcc_max_y & -bcc_max_z
c8_wfuel_reg = -c8_peb_in & c8_corner_bound
c8_nofuel_reg = -c8_peb_out & +c8_peb_in & c8_corner_bound

c8_centers = openmc.model.pack_spheres(triso_r[-1], region=-c8_peb_in, 
                                       num_spheres=19000, seed = 658469875)
c8_trisos = [openmc.model.TRISO(triso_r[-1], bg_triso_univ, center) 
             for center in c8_centers]

c8_wfuel = openmc.Cell(region=c8_wfuel_reg)
lower_left, upper_right = c8_wfuel.region.bounding_box
shape = (5, 5, 5)
pitch = (upper_right - lower_left)/shape
c8_lattice = openmc.model.create_triso_lattice(
    c8_trisos, lower_left, pitch, shape, graphite)

c8_wfuel.fill = c8_lattice
c8_nofuel = openmc.Cell(fill=graphite, region=c8_nofuel_reg)

c8_cells = [c8_wfuel, c8_nofuel]



bcc_out_bounds = ( +bcc_min_x & -bcc_max_x 
                  & +bcc_min_y & -bcc_max_y 
                  & +bcc_min_z & -bcc_max_z ) 
bcc_fill_reg = (bcc_out_bounds & +body_peb_out 
                & +c1_peb_out & +c2_peb_out 
                & +c3_peb_out & +c4_peb_out 
                & +c5_peb_out & +c6_peb_out 
                & +c7_peb_out & +c8_peb_out)
bcc_fill = openmc.Cell(region=bcc_fill_reg, fill=he)



cell_list = (body_cells + c1_cells + c2_cells + c3_cells + c4_cells 
             + c5_cells + c6_cells + c7_cells + c8_cells + [bcc_fill])
universe = openmc.Universe(cells=cell_list)

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

#materials = list(geometry.get_all_materials().values())
#openmc.Materials(materials).export_to_xml()

materials = openmc.Materials([ucoavg, ucot, buffer, pyc, sic, graphite, he])
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
#settings.run_mode = 'eigenvalue'
settings.verbosity = 6
settings.particles = 5000
settings.generations_per_batch = 5
settings.batches = 60
settings.inactive = 20
settings.seed = 987654321
settings.temperature = {'method' : 'interpolation', 'tolerance' : 10.0}
settings.output = {'tallies': False}
#settings.volume_calculations = [vol_calc]
settings.export_to_xml()

bcc_model = openmc.model.Model(geometry, materials, settings)

#chain file set as environment variable
uco_volume = 1.5209
ucot.volume = uco_volume

#Here is the coupled operator version
operator = openmc.deplete.CoupledOperator(bcc_model)

#but getting the micro xs and running with an independent operator might
#be faster:

#fluxes, micros = openmc.deplete.get_microxs_and_flux(bcc_model, materials)
#operator = openmc.deplete.IndependentOperator(materials, fluxes, micros)

# 1549 effective full power days over 6 passes = 6 258 day passes (~8.6 months)
# dt steps based on table 1 in adaptive burnup..., Walter and Manera
# 165/1549 gives that 1 EFPD is approximately 0.1 MWd/kgHM
d_steps = [1] + [4] + [4] + [10]*9 + [25]*10 + [50]*24
#d_steps = [1] #test to get to a good pcm

reactor_power = 165.0*(10**6) #165 MWth, converted to W
#220K pebs, 19k triso per peb, vol_kernel, uco density, wt percent of u in uco
uco_weight = (223000*19000*(4/3)*np.pi*triso_r[0]**3)*ucot.density*0.8945
specific_power = reactor_power/uco_weight #this is in W/gHM

celi = openmc.deplete.CELIIntegrator(

    operator, d_steps, power_density=specific_power, timestep_units = 'd')

celi.integrate()


