## Pebble Centroid Map Dataset v1.0.0

This dataset contains pebble centroid coordinates for the fuel pebbles in a
pebble-bed reactor.

Each dataset is paired with the input file that made it.  The list of dataset 
groups are as follows:

### Generic HTGR SMR 1
The "generic-start-pebs.txt" file is the original LAMMPS dump file with no
modifications.  It includes type 1 and type 2 pebbles, which are identical,
and were only used within LAMMPS to differentiate their color.  If only having 
one type of particle in LAMMPS is preferred, "generic-type1-pebs.txt" has been
modified to change all pebbles to type 1.  All particles are 3 cm (0.03 m) in
radius.  Note that the units of the coordinates and accompanying input file 
are in meters.

Files included:
- generic-start-pebs.txt
    - This file contains pebble centroid coordinates.  As this is a LAMMPS 
    dump file, it also includes a header containing information on the 
    timestep, number of pebbles, and the simulation's bounding box.
    Ready to use in a LAMMPS simulation with similar region and particle type
    definitions.
- generic-type1-pebs.txt
    - This file is identical to generic-start-pebs.txt, except that all
    pebbles in the file are type 1.
- generic-start-pebs.csv
    - This contains the same pebble centroid coordinates as the LAMMPS dump 
    file, but in csv file format and without the LAMMPS header, pebble ID, 
    or type columns.
- pour-220000-peb.txt
    - This file is the input file that created the pebble positions in LAMMPS.
- run_lammps.py
    - a simple script that runs pour-220000-peb.text.  It can also be run from
    the command line instead.  Additonally, the original run 
    used OMP_NUM_THREADS=4.
