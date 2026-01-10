import openmc
import openmc.deplete
import numpy as np
import random
#path_prefix = "~/ghastly/openmc_models/inf_lat_dep/"
dep_file = 'inf_lat_dep/dep-center/iter3/depletion_results.h5'
mat_file = 'inf_lat_dep/dep-center/iter3/materials.xml'
res = openmc.deplete.Results(dep_file)
dep_t = res.get_times()
step_comps = [res.export_to_materials(i, 
                                      path=mat_file)[0].get_nuclide_densities()
              for i in range(len(dep_t))]

#materials
#also double check isotopic compostions conventions in mats
#graphite based on a3-3, triso layers pulled from reported values in
#Neutronics characteristics of a 165 MWth Xe-100 reactor, Mulder et al


#for k-eigenvalue full core, we have 0-6 pass comps

#print(dep_t[([17, 20, 21, 22, 23, 24])])
#pass1 is from 1 to 17, inclusive. pass2 is 18, 19, and 20.  
#pass 3, 4, 5, and 6 are 21, 22, 23, and 24, respectively
#These are all going to change with the new, finer timesteps, but you *should*
#be able to predict them.  Also, now all passes will need to be over a range
#like pass 1 and 2 (but the bcc would use this so these should match)

#for now, use a set of dummy indices that you'll fix later
temp_comp = {}
for step in step_comps[1:4]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass1 = {}
for k, v in temp_comp.items():
    pass1[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[4:8]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1      
pass2 = {}
for k, v in temp_comp.items():
    pass2[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[8:12]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass3 = {}
for k, v in temp_comp.items():
    pass3[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[12:16]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass4 = {}
for k, v in temp_comp.items():
    pass4[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[16:20]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass5 = {}
for k, v in temp_comp.items():
    pass5[k] = v['iso']/v['count']

temp_comp ={}
for step in step_comps[20:]:
    for k, v in step.items():
        if k in temp_comp:
            temp_comp[k]['iso'] += v[1]
            temp_comp[k]['count'] += 1
        else:
            temp_comp[k] = {}
            temp_comp[k]['iso'] = v[1]
            temp_comp[k]['count'] = 1

pass6 = {}
for k, v in temp_comp.items():
    pass6[k] = v['iso']/v['count']


uco0 = openmc.Material(name='UCO_0')
uco0.set_density('g/cm3', 10.4)
uco0.add_nuclide("U235", 0.1386, percent_type='wo')
uco0.add_nuclide("U238",0.7559, percent_type='wo')
uco0.add_element("O", 0.06025, percent_type='wo')
uco0.add_element('C', 0.04523, percent_type='wo')
uco0.add_s_alpha_beta('c_Graphite')
#uco0.depletable = True
uco0.temperature = 1159.15 #K

uco1= openmc.Material(name='UCO_1')
uco1.set_density('g/cm3', 10.4)
uco1.add_components(pass1, percent_type = 'ao')
uco1.add_s_alpha_beta('c_Graphite')
#uco1.depletable = False
uco1.temperature = 1159.15 #K

uco2= openmc.Material(name='UCO_2')
uco2.set_density('g/cm3', 10.4)
uco2.add_components(pass2, percent_type = 'ao')
uco2.add_s_alpha_beta('c_Graphite')
#uco2.depletable = False
uco2.temperature = 1159.15 #K

uco3= openmc.Material(name='UCO_3')
uco3.set_density('g/cm3', 10.4)
uco3.add_components(pass3, percent_type = 'ao')
uco3.add_s_alpha_beta('c_Graphite')
#uco3.depletable = False
uco3.temperature = 1159.15 #K

uco4= openmc.Material(name='UCO_4')
uco4.set_density('g/cm3', 10.4)
uco4.add_components(pass4, percent_type = 'ao')
uco4.add_s_alpha_beta('c_Graphite')
#uco4.depletable = False
uco4.temperature = 1159.15 #K

uco5= openmc.Material(name='UCO_5')
uco5.set_density('g/cm3', 10.4)
uco5.add_components(pass5, percent_type = 'ao')
uco5.add_s_alpha_beta('c_Graphite')
#uco5.depletable = False
uco5.temperature = 1159.15 #K

uco6= openmc.Material(name='UCO_6')
uco6.set_density('g/cm3', 10.4)
uco6.add_components(pass6, percent_type = 'ao')
uco6.add_s_alpha_beta('c_Graphite')
#uco6.depletable = False
uco6.temperature = 1159.15 #K

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
sic.temperature = 1159.15 #K


graphite = openmc.Material(name='graphite')
graphite.set_density('kg/m3', 1700)
graphite.add_element('C', 1.0, percent_type='ao')
graphite.add_s_alpha_beta('c_Graphite')


he = openmc.Material(name='He')
he.set_density('atom/b-cm', 0.0006)
he.add_element('He', 1.0, percent_type='ao')
he.temperature = 778.15 #K

#from asm ss304 properties sheet
ss304 = openmc.Material(name='SS304')
ss304.set_density('g/cm3', 8.0)
ss304.add_element('Fe', 70.35, percent_type='wo')
ss304.add_element('Cr', 19.0, percent_type='wo')
ss304.add_element('Ni', 9.0, percent_type='wo')
ss304.add_element('Mn', 1.0, percent_type='wo')
ss304.add_element('Si', 0.5, percent_type='wo')
ss304.add_element('C', 0.05, percent_type='wo')
ss304.add_element('P', 0.04, percent_type='wo')
ss304.add_element('S', 0.03, percent_type='wo')

#geometry parameters - remember LAMMPS coords are in meters, this is in cm
#remember, triso order: kernel, buffer, pyc, sic, pyc
peb_or = 3.0 #outer radius of whole pebble
peb_ir = 2.5 # radius of the region that has trisos in it only
triso_r = [0.02125, 0.03125, 0.03525, 0.03875, 0.04275]

#active core params.  bottom of main cone(mcone) is at (0, 0, 0)
icyl_r = 24
icyl_h = 10
icone_h = 10
icone_t = 24
icone_b = 120
mcyl_r = 120
mcyl_h = 891
mcone_h = 54
mcone_t = 120
mcone_b = 24
mcone_zb = 0.0 #we'll orient things later based on the main cone exit
ocyl_r = 24
ocyl_h = 134

#reflector params
refl_side_t = 91
refl_top_t = 85
refl_bot_t = 36

#rcs params.  rod length = 6.6 m
rcs_r = 13/2

#he riser params
to_riser = 52 + (18/2) # distance from core center to riser center
riser_r = 9

#gap, barrel, and rpv params
gap1_t = 4 #between reflector and barrel
barrel_t = 4
gap2_t = 8 #between barrel and rpv
rpv_t = 9

#the triso bounds can stay the same, but the cell and univ definitions 
#need to be unique for every fuel comp

triso_bounds = [openmc.Sphere(r=r) for r in triso_r[:-1]]

triso_cells_0 = [openmc.Cell(fill=uco0, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_0 = openmc.Universe(cells=triso_cells_0)

triso_cells_1 = [openmc.Cell(fill=uco1, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_1 = openmc.Universe(cells=triso_cells_1)

triso_cells_2 = [openmc.Cell(fill=uco2, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_2 = openmc.Universe(cells=triso_cells_2)

triso_cells_3 = [openmc.Cell(fill=uco3, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_3 = openmc.Universe(cells=triso_cells_3)

triso_cells_4 = [openmc.Cell(fill=uco4, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_4 = openmc.Universe(cells=triso_cells_4)

triso_cells_5 = [openmc.Cell(fill=uco5, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_5 = openmc.Universe(cells=triso_cells_5)

triso_cells_6 = [openmc.Cell(fill=uco6, region=-triso_bounds[0]),
               openmc.Cell(fill=buffer,
                           region=+triso_bounds[0] & -triso_bounds[1]),
               openmc.Cell(fill=pyc, 
                           region=+triso_bounds[1] & -triso_bounds[2]),
               openmc.Cell(fill=sic, 
                           region=+triso_bounds[2] & -triso_bounds[3]),
               openmc.Cell(fill=pyc, region=+triso_bounds[3])]
triso_univ_6 = openmc.Universe(cells=triso_cells_6)

#you also need a unique triso pack/pebble inside definition for 0-6 pass, 
#but the sphere surfaces shouldn't need
#to be redefined each time
peb_in = openmc.Sphere(r = peb_ir)
peb_out = openmc.Sphere(r=peb_or)

# 0, aka fresh:
p0_wf_reg = -peb_in
p0_nof_reg = -peb_out & +peb_in

p0_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p0_wf_reg,
                                         num_spheres=19000, seed = 656574651)
p0_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_0, center) 
               for center in p0_triso_centers]

p0_wf = openmc.Cell(region=p0_wf_reg)
p0_nof = openmc.Cell(fill = graphite, region = p0_nof_reg)
lower_left, upper_right = p0_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p0_lattice = openmc.model.create_triso_lattice(
    p0_trisos, lower_left, pitch, shape, graphite)
p0_wf.fill = p0_lattice
p0_cells = [p0_wf, p0_nof]
p0_univ = openmc.Universe(cells = p0_cells)



#pass 1:
p1_wf_reg = -peb_in
p1_nof_reg = -peb_out & +peb_in

p1_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p1_wf_reg,
                                         num_spheres=19000, seed = 858932478)
p1_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_1, center) 
               for center in p1_triso_centers]

p1_wf = openmc.Cell(region=p1_wf_reg)
p1_nof = openmc.Cell(fill = graphite, region = p1_nof_reg)
lower_left, upper_right = p1_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p1_lattice = openmc.model.create_triso_lattice(
    p1_trisos, lower_left, pitch, shape, graphite)
p1_wf.fill = p1_lattice
p1_cells = [p1_wf, p1_nof]
p1_univ = openmc.Universe(cells = p1_cells)



#pass 2:
p2_wf_reg = -peb_in
p2_nof_reg = -peb_out & +peb_in

p2_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p2_wf_reg,
                                         num_spheres=19000, seed = 256485935)
p2_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_2, center) 
               for center in p2_triso_centers]

p2_wf = openmc.Cell(region=p2_wf_reg)
p2_nof = openmc.Cell(fill = graphite, region = p2_nof_reg)
lower_left, upper_right = p2_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p2_lattice = openmc.model.create_triso_lattice(
    p2_trisos, lower_left, pitch, shape, graphite)
p2_wf.fill = p2_lattice
p2_cells = [p2_wf, p2_nof]
p2_univ = openmc.Universe(cells = p2_cells)



#pass 3:

p3_wf_reg = -peb_in
p3_nof_reg = -peb_out & +peb_in

p3_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p3_wf_reg,
                                         num_spheres=19000, seed = 552136785)
p3_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_3, center) 
               for center in p3_triso_centers]

p3_wf = openmc.Cell(region=p3_wf_reg)
p3_nof = openmc.Cell(fill = graphite, region = p3_nof_reg)
lower_left, upper_right = p3_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p3_lattice = openmc.model.create_triso_lattice(
    p3_trisos, lower_left, pitch, shape, graphite)
p3_wf.fill = p3_lattice
p3_cells = [p3_wf, p3_nof]
p3_univ = openmc.Universe(cells = p3_cells)



#pass 4:

p4_wf_reg = -peb_in
p4_nof_reg = -peb_out & +peb_in

p4_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p4_wf_reg,
                                         num_spheres=19000, seed = 456375827)
p4_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_4, center) 
               for center in p4_triso_centers]

p4_wf = openmc.Cell(region=p4_wf_reg)
p4_nof = openmc.Cell(fill = graphite, region = p4_nof_reg)
lower_left, upper_right = p4_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p4_lattice = openmc.model.create_triso_lattice(
    p4_trisos, lower_left, pitch, shape, graphite)
p4_wf.fill = p4_lattice
p4_cells = [p4_wf, p4_nof]
p4_univ = openmc.Universe(cells = p4_cells)



#pass 5:

p5_wf_reg = -peb_in
p5_nof_reg = -peb_out & +peb_in

p5_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p5_wf_reg,
                                         num_spheres=19000, seed = 734844387)
p5_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_5, center) 
               for center in p5_triso_centers]

p5_wf = openmc.Cell(region=p5_wf_reg)
p5_nof = openmc.Cell(fill = graphite, region = p5_nof_reg)
lower_left, upper_right = p5_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p5_lattice = openmc.model.create_triso_lattice(
    p5_trisos, lower_left, pitch, shape, graphite)
p5_wf.fill = p5_lattice
p5_cells = [p5_wf, p5_nof]
p5_univ = openmc.Universe(cells = p5_cells)


#pass 6:

p6_wf_reg = -peb_in
p6_nof_reg = -peb_out & +peb_in

p6_triso_centers = openmc.model.pack_spheres(triso_r[-1], region=p6_wf_reg,
                                         num_spheres=19000, seed = 464865437)
p6_trisos = [openmc.model.TRISO(triso_r[-1], triso_univ_6, center) 
               for center in p6_triso_centers]

p6_wf = openmc.Cell(region=p6_wf_reg)
p6_nof = openmc.Cell(fill = graphite, region = p6_nof_reg)
lower_left, upper_right = p6_wf.region.bounding_box
shape = (3, 3, 3) 
pitch = (upper_right - lower_left)/shape
p6_lattice = openmc.model.create_triso_lattice(
    p6_trisos, lower_left, pitch, shape, graphite)
p6_wf.fill = p6_lattice
p6_cells = [p6_wf, p6_nof]
p6_univ = openmc.Universe(cells = p6_cells)

# now we need to read in the coords csv, and assign each coord a p#_univ.
# p0 can only be in the top half, p6 only in the bottom.

coord_file = "../examples/pre-made-lammps/generic-start-pebs.csv"
#lammps is in m, convert to cm for openmc
peb_coords = 100*np.loadtxt(coord_file, delimiter = ",", skiprows=1)
#now we need to split up the coords into their materials, based on the
#value of their z coordinate.  We want to have fresh pebbles only in the top
#half of the *cylindrical* portion of the main core.  so the cuttoff is:

halfway = mcone_zb + mcone_h + (mcyl_h/2)
#split the pebbles into top and bottom, then randomly assign compositions
#afterward
top_pebs = []
bot_pebs = []
for coord in peb_coords:
    if coord[2] > halfway:
        top_pebs += [(coord)]

    else:
        bot_pebs += [(coord)]
#awful name.  list where each entry is a dict with each key: value pair
#being the coord:univ pair for a specific pebble.
#NOTE TO FUTURE SELF:  WHEN YOU IMPLEMENT THIS AUTOMATICALLY, THIS INFO WILL
#NEED TO BE IN THE PEBBLE CLASS SOMEHOW, ALONGSIDE THE UNIQUE ID FROM LAMMPS
peb_params = []
top_univs = [p0_univ, p1_univ, p2_univ, p3_univ, p4_univ, p5_univ]
bot_univs = [p1_univ, p2_univ, p3_univ, p4_univ, p5_univ, p6_univ]
#LAMMPS organizes its particle coordinates.  This is helpful when moving them,
#but when we want to randomly assign materials, we can't step through the list
#or we might end up with an ordered array of pebbles we don't intend.
#soln: randomly shuffle a list of *index* positions to go over
num_bins = 6
bin_size = len(top_pebs)//num_bins
top_sample = random.sample(range(len(top_pebs)), len(top_pebs))
for i in range(num_bins):
    start = i*bin_size
    if i < num_bins-1:
        stop = (i+1)*bin_size
        for peb in top_sample[start:stop]:
            peb_params += [(top_univs[i], top_pebs[peb])]

    else:
        for peb in top_sample[start:]:
            peb_params += [(top_univs[i], top_pebs[peb])]

bot_sample = random.sample(range(len(bot_pebs)), len(bot_pebs))
for i in range(num_bins):
    start = i*bin_size
    if i < num_bins-1:
        stop = (i+1)*bin_size
        for peb in bot_sample[start:stop]:
            peb_params += [(bot_univs[i], bot_pebs[peb])]
    else:
        for peb in bot_sample[start:]:
            peb_params += [(bot_univs[i], bot_pebs[peb])]


#intake
icyl_side = openmc.ZCylinder(r = icyl_r)
icyl_top = openmc.ZPlane(z0 = (mcone_zb + mcone_h + mcyl_h + icone_h + icyl_h))
icyl_bot = openmc.ZPlane(z0 = (mcone_zb + mcone_h + mcyl_h + icone_h))
icone_H = -(icone_b*icone_h)/(icone_t-icone_b)
icone_z0 = mcone_zb + mcone_h + mcyl_h + icone_H
icone_r2 = (icone_b/icone_H)**2
icone_side = openmc.ZCone(z0 = icone_z0, r2 = icone_r2)

#main
mcyl_side = openmc.ZCylinder(r = mcyl_r)
mcyl_top = openmc.ZPlane(z0 = (mcone_zb + mcone_h + mcyl_h))
mcyl_bot = openmc.ZPlane(z0 = (mcone_zb + mcone_h))
#need to figure out where a full cone's apex would be, in order 
#to get the bottom where it's needed
mcone_H = -(mcone_t*mcone_h)/(mcone_b-mcone_t)
mcone_z0 = mcone_h - mcone_H
mcone_r2 = (mcone_t/mcone_H)**2
mcone_side = openmc.ZCone(z0 = mcone_z0, r2 = mcone_r2)
mcone_bot = openmc.ZPlane(z0 = mcone_zb)

#outtake
ocyl_side = openmc.ZCylinder(r = ocyl_r)
ocyl_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h))


icyl_reg = -icyl_side & -icyl_top & +icyl_bot
icone_reg = -icone_side & -icyl_bot & +mcyl_top
mcyl_reg = -mcyl_side & -mcyl_top & +mcyl_bot
mcone_reg = -mcone_side & -mcyl_bot & +mcone_bot
ocyl_reg = -ocyl_side & -mcone_bot & +ocyl_bot
active_core_reg = icyl_reg | icone_reg | mcyl_reg | mcone_reg | ocyl_reg

pebs = [openmc.model.TRISO(peb_or, peb[0], peb[1]) for peb in peb_params ]

active_core = openmc.Cell(region=active_core_reg)
#you can't automatically get the bounds of a region w/cones, but if it fits
#the cylinders, it should fit the cones, so:
bbox_reg = icyl_reg | mcyl_reg | ocyl_reg
bbox = openmc.Cell(region = bbox_reg)
lower_left_core, upper_right_core = bbox.region.bounding_box
shape_core = (7, 7, 7)
pitch_core = (upper_right_core - lower_left_core)/shape_core
core_lattice = openmc.model.create_triso_lattice(
    pebs, lower_left_core, pitch_core, shape_core, he)

active_core.fill=core_lattice

refl_side = openmc.ZCylinder(r = (mcyl_r + refl_side_t))
refl_top = openmc.ZPlane(z0 = (mcone_zb + mcone_h + mcyl_h + 
                               icone_h + icyl_h + refl_top_t))
refl_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h - refl_bot_t))

#add in 24 helium risers, 18 RCS (9 normal, 9 SCRAM)
#he first.  from SAM model, the close edge of the riser is 52 cm into the
#reflector, then the riser diameter is 18 cm, leaving 21 cm between the
#far edge of the riser and the outer graphite reflector boundary

he_rad = (np.pi/12)*np.linspace(0,23,24)

he_side = [openmc.ZCylinder(x0 = to_riser*np.cos(rad),
                            y0 = to_riser*np.sin(rad),
                            r = riser_r) for rad in he_rad]
he_reg = (-he_side[0] & 
          -he_side[1] &
          -he_side[2] &
          -he_side[3] &
          -he_side[4] &
          -he_side[5] &
          -he_side[6] &
          -he_side[7] &
          -he_side[8] &
          -he_side[9] &
          -he_side[10] &
          -he_side[11] &
          -he_side[12] &
          -he_side[13] &
          -he_side[14] &
          -he_side[15] &
          -he_side[16] &
          -he_side[17] &
          -he_side[18] &
          -he_side[19] &
          -he_side[20] &
          -he_side[21] &
          -he_side[22] &
          -he_side[23] &
          -icyl_top & +ocyl_bot)
he_risers = openmc.Cell(region = he_reg, fill = he)


refl_reg = -refl_side & -refl_top & + refl_bot & ~active_core_reg & ~he_reg
refl = openmc.Cell(region = refl_reg, fill = graphite)

#gap1
gap1_side = openmc.ZCylinder(r = mcyl_r + refl_side_t + gap1_t)
gap1_top = openmc.ZPlane(z0 =(mcone_zb + mcone_h + mcyl_h +
                              icone_h + icyl_h + refl_top_t + gap1_t))
gap1_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h - refl_bot_t - gap1_t))
gap1_side_reg = -gap1_side & +refl_side & -refl_top & +refl_bot
gap1_top_reg = -gap1_side & -gap1_top & +refl_top
gap1_bot_reg = -gap1_side & -refl_bot & +gap1_bot
gap1_reg = gap1_side_reg | gap1_top_reg | gap1_bot_reg
gap1 = openmc.Cell(region = gap1_reg, fill = he)

#barrel
barrel_side = openmc.ZCylinder(r = (mcyl_r + refl_side_t + 
                                    gap1_t + barrel_t))
barrel_top = openmc.ZPlane(z0 =(mcone_zb + mcone_h + mcyl_h +
                              icone_h + icyl_h + refl_top_t + 
                                gap1_t + barrel_t))
barrel_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h - refl_bot_t - 
                                 gap1_t - barrel_t))
barrel_side_reg = -barrel_side & +gap1_side & -refl_top & +refl_bot
barrel_top_reg = -barrel_side & -barrel_top & +refl_top & ~gap1_top_reg
barrel_bot_reg = -barrel_side & -refl_bot & +barrel_bot & ~gap1_bot_reg
barrel_reg = barrel_side_reg | barrel_top_reg | barrel_bot_reg
barrel = openmc.Cell(region = barrel_reg, fill = ss304)

#gap2
gap2_side = openmc.ZCylinder(r = (mcyl_r + refl_side_t + 
                                    gap1_t + barrel_t + gap2_t))
gap2_top = openmc.ZPlane(z0 =(mcone_zb + mcone_h + mcyl_h +
                              icone_h + icyl_h + refl_top_t + 
                                gap1_t + barrel_t + gap2_t))
gap2_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h - refl_bot_t - 
                                 gap1_t - barrel_t - gap2_t))
gap2_side_reg = -gap2_side & +barrel_side & -refl_top & +refl_bot
gap2_top_reg = (-gap2_side & -gap2_top & +refl_top & 
                ~gap1_top_reg & ~barrel_top_reg)
gap2_bot_reg = (-gap2_side & -refl_bot & +gap2_bot &
                ~gap1_bot_reg & ~barrel_bot_reg)
gap2_reg = gap2_side_reg | gap2_top_reg | gap2_bot_reg
gap2 = openmc.Cell(region = gap2_reg, fill = he)

#rpv
rpv_side = openmc.ZCylinder(r = (mcyl_r + refl_side_t + 
                                 gap1_t + barrel_t + gap2_t + rpv_t),
                            boundary_type='vacuum')
rpv_top = openmc.ZPlane(z0 =(mcone_zb + mcone_h + mcyl_h +
                              icone_h + icyl_h + refl_top_t + 
                                gap1_t + barrel_t + gap2_t + rpv_t),
                        boundary_type='vacuum')
rpv_bot = openmc.ZPlane(z0 = (mcone_zb - ocyl_h - refl_bot_t - 
                                 gap1_t - barrel_t - gap2_t - rpv_t),
                        boundary_type='vacuum')
rpv_side_reg = -rpv_side & +gap2_side & -refl_top & +refl_bot
rpv_top_reg = (-rpv_side & -rpv_top & +refl_top & 
                ~gap1_top_reg & ~barrel_top_reg & ~gap2_top_reg)
rpv_bot_reg = (-rpv_side & -refl_bot & +rpv_bot &
                ~gap1_bot_reg & ~barrel_bot_reg & ~gap2_bot_reg)
rpv_reg = rpv_side_reg | rpv_top_reg | rpv_bot_reg
rpv = openmc.Cell(region = rpv_reg, fill = ss304)



universe = openmc.Universe(cells=[active_core, refl, gap1, barrel, gap2, rpv])

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

materials = list(geometry.get_all_materials().values())
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.run_mode = 'plot'
settings.export_to_xml()


xz_cell=openmc.Plot()
xz_cell.basis='xz'
xz_cell.origin=(0,0, (-ocyl_h - refl_bot_t + mcone_h + 
                      mcyl_h + icone_h + icyl_h+ refl_top_t)/2 )
xz_cell.width = (440, 1250)
xz_cell.pixels = (880, 2500)
xz_cell.color_by='cell'


xz_mats=openmc.Plot()
xz_mats.basis='xz'
xz_mats.origin=(0,0, (-ocyl_h - refl_bot_t + mcone_h + 
                      mcyl_h + icone_h + icyl_h+ refl_top_t)/2 )
xz_mats.width = (440, 1250)
xz_mats.pixels = (880, 2500)
xz_mats.color_by='material'
xz_mats.colors = {graphite: 'gray', he: 'lightblue', buffer: 'mediumslateblue', 
                  pyc: 'lavender', sic: 'darkslateblue',
                  uco0: 'lavenderblush', uco1:'mistyrose' , 
                  uco2:'pink' , uco3:'lightpink' , 
                  uco4:'palevioletred' , uco5:'crimson' , uco6:'maroon'}



plots = openmc.Plots([xz_cell, xz_mats])
plots.export_to_xml()
openmc.plot_geometry()
