import openmc
import openmc.deplete
import numpy as np

#This script generates .xml files that all model iterations will use.

peb_or = 3.0 #outer radius of whole pebble
peb_ir = 2.5 # radius of the region that has trisos in it only
triso_r = [0.02125, 0.03125, 0.03525, 0.03875, 0.04275]
bcc_l = 100/(2697**(1/3)) #based on pf = 5394 pebs/m3
c_coord = 0.5*bcc_l

#materials
#graphite based on a3-3, triso layers pulled from reported values in
#Neutronics characteristics of a 165 MWth Xe-100 reactor, Mulder et al

ucodep = openmc.Material(name='UCO_DEP', material_id=1)
ucodep.set_density('g/cm3', 10.4)
ucodep.add_nuclide("U235", 0.1386, percent_type='wo')
ucodep.add_nuclide("U238", 0.7559, percent_type='wo')
ucodep.add_element("O", 0.06025, percent_type='wo')
ucodep.add_element('C', 0.04523, percent_type='wo')
ucodep.add_s_alpha_beta('c_Graphite')
ucodep.depletable = True
ucodep.temperature = 1159.15 #K
ucodep.volume = 19000*(4/3)*np.pi*triso_r[0]**3

#placeholder for defining geometry
ucoavg = openmc.Material(name='UCO_AVG', material_id=2)

buffer = openmc.Material(name='buffer', material_id=3)
buffer.set_density('g/cm3', 1.05)
buffer.add_element('C', 1.0, percent_type='ao')
buffer.add_s_alpha_beta('c_Graphite')
buffer.temperature = 1159.15 #K

pyc = openmc.Material(name='PyC', material_id=4)
pyc.set_density('g/cm3', 1.9)
pyc.add_element('C', 1.0, percent_type='ao')
pyc.add_s_alpha_beta('c_Graphite')
pyc.temperature = 1159.15 #K

sic = openmc.Material(name='SiC', material_id=5)
sic.set_density('g/cm3', 3.2)
sic.add_element('C', 0.5, percent_type='ao')
sic.add_element('Si', 0.5, percent_type='ao')
sic.add_s_alpha_beta('c_Graphite')

graphite = openmc.Material(name='graphite', material_id=6)
graphite.set_density('kg/m3', 1700)
graphite.add_element('C', 1.0, percent_type='ao')
graphite.add_s_alpha_beta('c_Graphite')

he = openmc.Material(name='He', material_id=7)
he.set_density('atom/b-cm', 0.0006)
he.add_element('He', 1.0, percent_type='ao')
he.temperature = 778.15 #K


#define the edges of the unit cell
bcc_min_x = openmc.XPlane(x0=-c_coord, boundary_type='reflective')
bcc_max_x = openmc.XPlane(x0=c_coord, boundary_type='reflective')
bcc_min_y = openmc.YPlane(y0=-c_coord, boundary_type='reflective')
bcc_max_y = openmc.YPlane(y0=c_coord, boundary_type='reflective')
bcc_min_z = openmc.ZPlane(z0=-c_coord, boundary_type='reflective')
bcc_max_z = openmc.ZPlane(z0=c_coord, boundary_type='reflective')

#define the cells and universe for the triso particles
triso_reg = [openmc.Sphere(r=r) for r in triso_r[:-1]]
tr_triso_cells = [openmc.Cell(fill=ucodep, region=-triso_reg[0]),
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

seeds = [978397880, 987432789, 895490889, 
         356657913, 353623684, 897578459, 
         429875621, 489638795, 658469875]


#center pebble
body_peb_in = openmc.Sphere(r = peb_ir)
body_wfuel_bound = -body_peb_in
body_peb_out = openmc.Sphere(r=peb_or)
body_nofuel_reg = -body_peb_out & +body_peb_in

body_centers = openmc.model.pack_spheres(triso_r[-1], region=body_wfuel_bound,
                                         num_spheres=19000, seed = seeds[0])
body_trisos = [openmc.model.TRISO(triso_r[-1], tr_triso_univ, center) 
               for center in body_centers]

body_wfuel = openmc.Cell(region=body_wfuel_bound)
lower_left, upper_right = body_wfuel.region.bounding_box
shape = (10, 10, 10)
pitch = (upper_right - lower_left)/shape
body_lattice = openmc.model.create_triso_lattice(
    body_trisos, lower_left, pitch, shape, graphite)

body_wfuel.fill = body_lattice
body_nofuel = openmc.Cell(fill=graphite, region=body_nofuel_reg)

body_cells = [body_wfuel, body_nofuel]


# -1,-1,-1 corner
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
                                       num_spheres=19000, seed = seeds[1])
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
                                       num_spheres=19000, seed = seeds[2])
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
                                       num_spheres=19000, seed = seeds[3])
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
                                       num_spheres=19000, seed = seeds[4])
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
                                       num_spheres=19000, seed = seeds[5])
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
                                       num_spheres=19000, seed = seeds[6])
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
                                       num_spheres=19000, seed = seeds[7])
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
                                       num_spheres=19000, seed = seeds[8])
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

materials = openmc.Materials([ucodep, buffer, pyc, sic, graphite, he])
openmc.Materials(materials).export_to_xml()

settings = openmc.Settings()
settings.verbosity = 7
settings.particles = 5000
settings.generations_per_batch = 10
settings.batches = 75
settings.inactive = 25
settings.seed = 987654321
settings.temperature = {'method' : 'interpolation', 'tolerance' : 10.0}
settings.output = {'tallies': False,
                   'summary' : False}
settings.export_to_xml()


