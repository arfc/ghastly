# Ghastly

Ghastly is a python-based package to externally couple a Discrete Element Method (DEM) simulation program - to model pebble bed flow patterns - to a Monte Carlo code for reactor physics and depletion simulations.

Ghastly input files are in JSON format, and contain the various parameters both [LAMMPS](https://www.lammps.org/) and [OpenMC](https://openmc.org/) need to define physics parameters and active core geometries.  Ghastly interacts with LAMMPS purely by creating input files for LAMMPS and using data from LAMMPS output files to inform analysis. For OpenMC, Ghastly will automatically generate the material and geometry sections for the active core - but users will need to supply their own "shell" to place the active core in, usually encompassing everything from the reflector outward. The examples directory has simple input files, including models of a generic High Temperature Gas-Cooled Reactor (HTGR) Small Modular Reactor (SMR). The default jinja templates are included in the ghastly/templates directory.

Each LAMMPS template has some basic settings hard-coded - these are generally options that either cannot be changed for Ghastly to work, options that are unlikely to need to be changed, or simply basic settings that shouldn't impact the results, such as particle colors. While Ghastly generally defaults to the pre-made templates, users ccan also create their own custom templates to suit their needs. Users unfamiliar with LAMMPS should read the documentation  carefully before making changes to any fix or command. 

NOTE: Ghastly assumes LAMMPS is using the SI definitions for its units. This means that distance in LAMMPS is in meters. However, OpenMC strictly uses centimeters. At this time, Ghastly automatically converts values from LAMMPS from meters to centimeters - so using centimeters in LAMMPS will cause the OpenMC dimensions to be off. The Ghastly input file should be in meters.

Ghastly is not packaged with either LAMMPS or OpenMC source code - it is recommended that users create a conda environment to contain LAMMPS, OpenMC, and Ghastly installations.  Installation advice/instructions and minimum requirements are provided below.

## Ghastly can:
- Fill both upwards and downwards flowing PBRs with an initial bed configuration
- Recirculate downwards flowing (HTGR) pebble beds.

## Ghastly cannot:
- Model PBRs with a central reflector
- Model PBRs with multiple inlet and outlet chutes

Currently, Ghastly supports coupling LAMMMPS to OpenMC only, but developing and submitting an update to add further code support via PR is welcome.

## Installation:

The included setup_script.sh has the basic commands to install OpenMC, LAMMPS, and Ghastly, along with their dependencies, inside a pre-existing conda-environment, but it can also be used as a guide for those that want to customize their install, or add it to an environment that has other software.

Note that the steps outlined in the setup script will install as many things as possible inside the conda environment - including software built from source, using ${CONDA_PREFIX}.  If your current setup already has many of the dependencies installed in the base environment, your system may not be able to properly find the necessary libraries or files inside the conda environment.

## LAMMPS Troubleshooting:
I would also recommend looking at LAMMPS's own troubleshooting page - for this sort of problem, here are some common issues, and potential solutions:

- Pebbles overlapping/'bursting' and escaping the simulation box:
    - When a particle moves too far before LAMMPS recalculates forces from external bodies such as neighbor pebbles and walls on it, it can overlap, causing the particle to escape the simulation bounds, or have a repulsive force so great it 'rockets' away from the neighbor pebble
        - The timestep in Ghastly is set, by default, to be 1/50 of the collision time, based on guidance from the LAMMPS documentation page and [3 papers LAMMPS mentions].  It's generally recommended that the timestep remain at this value - and anything over a microsecond is likely too large.
        - The skin value is the distance a particle must move before LAMMPS recalculates forces acting on it.  If this is too high, a particle may overlap a neightbor before the simulation catches it.
        - Ghastly includes a line to explicitly set the pebble size - if you are using a custom template that does not do this, it's possible that LAMMPS is using the default particle size (which will be far too large for a PBR, in SI units)
- Pebbles get stuck between regions:
    - Check that the top and bottom of the regions you've defined are open.  For a downwards flowing system, only the very bottom of the outlet region should be closed.  For an upwards flowing region this is still the case - but remember that the outlet region should be at the geometric 'top' of the model.

