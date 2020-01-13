Input script structure
======================

This page describes the structure of a typical LAMMPS input script.
The examples directory in the LAMMPS distribution contains many sample
input scripts; it is discussed on the :doc:`Examples <Examples>` doc
page.

A LAMMPS input script typically has 4 parts:

1. Initialization
2. Atom definition
3. Settings
4. Run a simulation

The last 2 parts can be repeated as many times as desired.  I.e. run a
simulation, change some settings, run some more, etc.  Each of the 4
parts is now described in more detail.  Remember that almost all
commands need only be used if a non-default value is desired.

(1) Initialization

Set parameters that need to be defined before atoms are created or
read-in from a file.

The relevant commands are :doc:`units <units>`,
:doc:`dimension <dimension>`, :doc:`newton <newton>`,
:doc:`processors <processors>`, :doc:`boundary <boundary>`,
:doc:`atom\_style <atom_style>`, :doc:`atom\_modify <atom_modify>`.

If force-field parameters appear in the files that will be read, these
commands tell LAMMPS what kinds of force fields are being used:
:doc:`pair\_style <pair_style>`, :doc:`bond\_style <bond_style>`,
:doc:`angle\_style <angle_style>`, :doc:`dihedral\_style <dihedral_style>`,
:doc:`improper\_style <improper_style>`.

(2) Atom definition

There are 3 ways to define atoms in LAMMPS.  Read them in from a data
or restart file via the :doc:`read\_data <read_data>` or
:doc:`read\_restart <read_restart>` commands.  These files can contain
molecular topology information.  Or create atoms on a lattice (with no
molecular topology), using these commands: :doc:`lattice <lattice>`,
:doc:`region <region>`, :doc:`create\_box <create_box>`,
:doc:`create\_atoms <create_atoms>`.  The entire set of atoms can be
duplicated to make a larger simulation using the
:doc:`replicate <replicate>` command.

(3) Settings

Once atoms and molecular topology are defined, a variety of settings
can be specified: force field coefficients, simulation parameters,
output options, etc.

Force field coefficients are set by these commands (they can also be
set in the read-in files): :doc:`pair\_coeff <pair_coeff>`,
:doc:`bond\_coeff <bond_coeff>`, :doc:`angle\_coeff <angle_coeff>`,
:doc:`dihedral\_coeff <dihedral_coeff>`,
:doc:`improper\_coeff <improper_coeff>`,
:doc:`kspace\_style <kspace_style>`, :doc:`dielectric <dielectric>`,
:doc:`special\_bonds <special_bonds>`.

Various simulation parameters are set by these commands:
:doc:`neighbor <neighbor>`, :doc:`neigh\_modify <neigh_modify>`,
:doc:`group <group>`, :doc:`timestep <timestep>`,
:doc:`reset\_timestep <reset_timestep>`, :doc:`run\_style <run_style>`,
:doc:`min\_style <min_style>`, :doc:`min\_modify <min_modify>`.

Fixes impose a variety of boundary conditions, time integration, and
diagnostic options.  The :doc:`fix <fix>` command comes in many flavors.

Various computations can be specified for execution during a
simulation using the :doc:`compute <compute>`,
:doc:`compute\_modify <compute_modify>`, and :doc:`variable <variable>`
commands.

Output options are set by the :doc:`thermo <thermo>`, :doc:`dump <dump>`,
and :doc:`restart <restart>` commands.

(4) Run a simulation

A molecular dynamics simulation is run using the :doc:`run <run>`
command.  Energy minimization (molecular statics) is performed using
the :doc:`minimize <minimize>` command.  A parallel tempering
(replica-exchange) simulation can be run using the
:doc:`temper <temper>` command.
