Input script structure
======================

This page describes the structure of a typical LAMMPS input script.
The examples directory in the LAMMPS distribution contains many sample
input scripts; it is discussed on the :doc:`Examples <Examples>` doc
page.

A LAMMPS input script typically has 4 parts:

1. :ref:`Initialization <init>`
2. :ref:`System definition <system>` 
3. :ref:`Simulation settings <settings>`
4. :ref:`Run a simulation <run>`

The last 2 parts can be repeated as many times as desired.  I.e. run a
simulation, change some settings, run some more, etc.  Each of the 4
parts is now described in more detail.  Remember that almost all
commands need only be used if a non-default value is desired.

.. _init:

Initialization
------------------------------

Set parameters that need to be defined before atoms are created or
read-in from a file.

The relevant commands are :doc:`units <units>`,
:doc:`dimension <dimension>`, :doc:`newton <newton>`,
:doc:`processors <processors>`, :doc:`boundary <boundary>`,
:doc:`atom_style <atom_style>`, :doc:`atom_modify <atom_modify>`.

If force-field parameters appear in the files that will be read, these
commands tell LAMMPS what kinds of force fields are being used:
:doc:`pair_style <pair_style>`, :doc:`bond_style <bond_style>`,
:doc:`angle_style <angle_style>`, :doc:`dihedral_style <dihedral_style>`,
:doc:`improper_style <improper_style>`.

.. _system:

System definition
------------------------------

There are 3 ways to define the simulation cell and reserve space for
force field info and fill it with atoms in LAMMPS.  Read them in from
(1) a data file or (2) a restart file via the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands,
respectively.  These files can also contain molecular topology
information.  Or (3) create a simulation cell and fill it with atoms on
a lattice (with no molecular topology), using these commands:
:doc:`lattice <lattice>`, :doc:`region <region>`, :doc:`create_box
<create_box>`, :doc:`create_atoms <create_atoms>` or
:doc:`read_dump <read_dump>`.

The entire set of atoms can be duplicated to make a larger simulation
using the :doc:`replicate <replicate>` command.

.. _settings:

Simulation settings
------------------------------

Once atoms and molecular topology are defined, a variety of settings
can be specified: force field coefficients, simulation parameters,
output options, and more.

Force field coefficients are set by these commands (they can also be
set in the read-in files): :doc:`pair_coeff <pair_coeff>`,
:doc:`bond_coeff <bond_coeff>`, :doc:`angle_coeff <angle_coeff>`,
:doc:`dihedral_coeff <dihedral_coeff>`,
:doc:`improper_coeff <improper_coeff>`,
:doc:`kspace_style <kspace_style>`, :doc:`dielectric <dielectric>`,
:doc:`special_bonds <special_bonds>`.

Various simulation parameters are set by these commands:
:doc:`neighbor <neighbor>`, :doc:`neigh_modify <neigh_modify>`,
:doc:`group <group>`, :doc:`timestep <timestep>`,
:doc:`reset_timestep <reset_timestep>`, :doc:`run_style <run_style>`,
:doc:`min_style <min_style>`, :doc:`min_modify <min_modify>`.

Fixes impose a variety of boundary conditions, time integration, and
diagnostic options.  The :doc:`fix <fix>` command comes in many flavors.

Various computations can be specified for execution during a
simulation using the :doc:`compute <compute>`,
:doc:`compute_modify <compute_modify>`, and :doc:`variable <variable>`
commands.

Output options are set by the :doc:`thermo <thermo>`, :doc:`dump <dump>`,
and :doc:`restart <restart>` commands.

.. _run:

Run a simulation
------------------------------

A molecular dynamics simulation is run using the :doc:`run <run>`
command.  Energy minimization (molecular statics) is performed using
the :doc:`minimize <minimize>` command.  A parallel tempering
(replica-exchange) simulation can be run using the
:doc:`temper <temper>` command.
