## Pebble Centroid Map Dataset v1.0.0

This dataset contains pebble centroid coordinates for the fuel pebbles in a
pebble-bed reactor.

Each dataset is paired with the input file that made it.  The list of dataset groups are as follows:

### Generic HTGR SMR 1
All particles, both type 1 and type 2, are 3[cm] or 0.03 [m] in radius.  Note that the units of the coordinates and accompanying input file are in meters.
Files included:
- generic-start-pebs.txt
    - contains pebble centroid coordinates.  As this is a LAMMPS dump file, it also includes a header containing information on the timestep, number of pebbles, and the simulation's bounding box.  Ready to use in a LAMMPS simulation with similar region and particle type definitions.
- generic-start-pebs.csv
    - contains the same pebble centroid coordinates as the LAMMPS dumpfile, but in csv file format and without the LAMMPS header, pebble ID, or type columns.
- pour-220000-peb.txt
    - input file that created the pebble positions in LAMMPS.
