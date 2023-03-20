.. index:: fix mdi/qm

fix mdi/qm command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mdi/qm keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* mdi/qm = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *virial* or *add* or *every* or *connect* or *elements*

  .. parsed-literal::

       *virial* args = *yes* or *no*
         yes = request virial tensor from server code
         no = do not request virial tensor from server code
       *add* args = *yes* or *no*
         yes = add returned value from server code to LAMMPS quantities
         no = do not add returned values to LAMMPS quantities
       *every* args = Nevery
         Nevery = request values from server code once every Nevery steps
       *connect* args = *yes* or *no*
         yes = perform a one-time connection to the MDI engine code
         no = do not perform the connection operation
       *elements* args = N_1 N_2 ... N_ntypes
         N_1,N_2,...N_ntypes = atomic number for each of ntypes LAMMPS atom types

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mdi/qm
   fix 1 all mdi/qm virial yes
   fix 1 all mdi/qm add no every 100 elements 13 29

Description
"""""""""""

.. versionadded:: 3Aug2022

This command enables LAMMPS to act as a client with another server
code that will compute the total energy, per-atom forces, and total
virial for atom conformations and simulation box size/shapes that
LAMMPS sends it.

Typically the server code will be a quantum mechanics (QM) code, hence
the name of the fix.  However this is not required, the server code
could be another classical molecular dynamics code or LAMMPS itself.
The server code must support use of the `MDI Library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_ as
explained below.

These are example use cases for this fix, discussed further below:

* perform an ab initio MD (AIMD) simulation with quantum forces
* perform an energy minimization with quantum forces
* perform a nudged elastic band (NEB) calculation with quantum forces
* perform a QM calculation for a series of independent systems which LAMMPS reads or generates

The code coupling performed by this command is done via the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_.
LAMMPS runs as an MDI driver (client), and sends MDI commands to an
external MDI engine code (server), e.g. a QM code which has support
for MDI.  See the :doc:`Howto mdi <Howto_mdi>` page for more
information about how LAMMPS can operate as either an MDI driver or
engine.

The examples/mdi directory contains input scripts using this fix in
the various use cases discussed below.  In each case, two instances of
LAMMPS are used, once as an MDI driver, once as an MDI engine
(surrogate for a QM code).  The examples/mdi/README file explains how
to launch two codes so that they communicate via the MDI library using
either MPI or sockets.  Any QM code that supports MDI could be used in
place of LAMMPS acting as a QM surrogate.  See the :doc:`Howto mdi
<Howto_mdi>` page for a current list (March 2022) of such QM codes.

Note that an engine code can support MDI in either or both of two
modes.  It can be used as a stand-alone code, launched at the same
time as LAMMPS.  Or it can be used as a plugin library, which LAMMPS
loads.  See the :doc:`mdi plugin <mdi>` command for how to trigger
LAMMPS to load a plugin library.  The examples/mdi/README file
explains how to launch the two codes in either mode.

----------

The *virial* keyword setting of yes or no determines whether
LAMMPS will request the QM code to also compute and return
a 6-element symmetric virial tensor for the system.

The *add* keyword setting of *yes* or *no* determines whether the
energy and forces and virial returned by the QM code will be added to
the LAMMPS internal energy and forces and virial or not.  If the
setting is *no* then the default :doc:`fix_modify energy <fix_modify>`
and :doc:`fix_modify virial <fix_modify>` settings are also set to
*no* and your input scripts should not set them to yes.  See more
details on these fix_modify settings below.

Whatever the setting for the *add* keyword, the QM energy, forces, and
virial will be stored by the fix, so they can be accessed by other
commands.  See details below.

The *every* keyword determines how often the QM code will be invoked
during a dynamics run with the current LAMMPS simulation box and
configuration of atoms.  The QM code will be called once every
*Nevery* timesteps.

The *connect* keyword determines whether this fix performs a one-time
connection to the QM code.  The default is *yes*.  The only time a
*no* is needed is if this command is used multiple times in an input
script.  E.g. if it used inside a loop which also uses the :doc:`clear
<clear>` command to destroy the system (including any defined fixes).
See the examples/mdi/in.series.driver script as an example of this,
where LAMMPS is using the QM code to compute energy and forces for a
series of system configurations.  In this use case *connect no*
is used along with the :doc:`mdi connect and exit <mdi>` command
to one-time initiate/terminate the connection outside the loop.

The *elements* keyword allows specification of what element each
LAMMPS atom type corresponds to.  This is specified by the atomic
number of the element, e.g. 13 for Al.  An atomic number must be
specified for each of the ntypes LAMMPS atom types.  Ntypes is
typically specified via the create_box command or in the data file
read by the read_data command.

If this keyword is specified, then this fix will send the MDI
">ELEMENTS" command to the engine, to ensure the two codes are
consistent in their definition of atomic species.  If this keyword is
not specified, then this fix will send the MDI >TYPES command to the
engine.  This is fine if both the LAMMPS driver and the MDI engine are
initialized so that the atom type values are consistent in both codes.

----------

The following 3 example use cases are illustrated in the examples/mdi
directory.  See its README file for more details.

(1) To run an ab initio MD (AIMD) dynamics simulation, or an energy
minimization with QM forces, or a multi-replica NEB calculation, use
*add yes* and *every 1* (the defaults).  This is so that every time
LAMMPS needs energy and forces, the QM code will be invoked.

Both LAMMPS and the QM code should define the same system (simulation
box, atoms and their types) in their respective input scripts.  Note
that on this scenario, it may not be necessary for LAMMPS to define a
pair style or use a neighbor list.

LAMMPS will then perform the timestepping or minimization iterations
for the simulation.  At the point in each timestep or iteration when
LAMMPS needs the force on each atom, it communicates with the engine
code.  It sends the current simulation box size and shape (if they
change dynamically, e.g. during an NPT simulation), and the current
atom coordinates.  The engine code computes quantum forces on each
atom and the total energy of the system and returns them to LAMMPS.

Note that if the AIMD simulation is an NPT or NPH model, or the energy
minimization includes :doc:`fix box relax <fix_box_relax>` to
equilibrate the box size/shape, then LAMMPS computes a pressure.  This
means the *virial* keyword should be set to *yes* so that the QM
contribution to the pressure can be included.

(2) To run dynamics with a LAMMPS interatomic potential, and evaluate
the QM energy and forces once every 1000 steps, use *add no* and
*every 1000*.  This could be useful for using an MD run to generate
randomized configurations which are then passed to the QM code to
produce training data for a machine learning potential.  A :doc:`dump
custom <dump>` command could be invoked every 1000 steps to dump the
atom coordinates and QM forces to a file.  Likewise the QM energy and
virial could be output with the :doc:`thermo_style custom
<thermo_style>` command.

(3) To do a QM evaluation of energy and forces for a series of *N*
independent systems (simulation box and atoms), use *add no* and
*every 1*.  Write a LAMMPS input script which loops over the *N*
systems.  See the :doc:`Howto multiple <Howto_multiple>` doc page for
details on looping and removing old systems.  The series of systems
could be initialized by reading them from data files with
:doc:`read_data <read_data>` commands.  Or, for example, by using the
:doc:`lattice <lattice>` , :doc:`create_atoms <create_atoms>`,
:doc:`delete_atoms <delete_atoms>`, and/or :doc:`displace_atoms
random <displace_atoms>` commands to generate a series of different
systems.  At the end of the loop perform :doc:`run 0 <run>` and
:doc:`write_dump <write_dump>` commands to invoke the QM code and
output the QM energy and forces.  As in (2) this be useful to produce
QM data for training a machine learning potential.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the potential energy computed by the QM code to the
global potential energy of the system as part of :doc:`thermodynamic
output <thermo_style>`.  The default setting for this fix is
:doc:`fix_modify energy yes <fix_modify>`, unless the *add* keyword is
set to *no*, in which case the default setting is *no*.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by
this fix to add the contribution computed by the QM code to the global
pressure of the system as part of :doc:`thermodynamic output
<thermo_style>`.  The default setting for this fix is :doc:`fix_modify
virial yes <fix_modify>`, unless the *add* keyword is set to *no*, in
which case the default setting is *no*.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the energy
returned by the QM code.  The scalar value calculated by this fix is
"extensive".

This fix also computes a global vector with of length 6 which contains
the symmetric virial tensor values returned by the QM code.  It can
likewise be accessed by various :doc:`output commands <Howto_output>`.

The ordering of values in the symmetric virial tensor is as follows:
vxx, vyy, vzz, vxy, vxz, vyz.  The values will be in pressure
:doc:`units <units>`.

This fix also computes a peratom array with 3 columns which contains
the peratom forces returned by the QM code.  It can likewise be
accessed by various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Assuming the *add* keyword is set to *yes* (the default), the forces
computed by the QM code are used during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the potential energy associated with the QM forces to
   be included in the total potential energy of the system (the
   quantity being minimized), you MUST not disable the
   :doc:`fix_modify <fix_modify>` *energy* option for this fix, which
   means the *add* keyword should also be set to *yes* (the default).


Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The QM code does not currently compute and return per-atom energy or
per-atom virial contributions.  So they will not show up as part of
the calculations performed by the :doc:`compute pe/atom
<compute_pe_atom>` or :doc:`compute stress/atom <compute_stress_atom>`
commands.

To use LAMMPS as an MDI driver in conjunction with other MDI-enabled
codes (MD or QM codes), the :doc:`units <units>` command should be
used to specify *real* or *metal* units.  This will ensure the correct
unit conversions between LAMMPS and MDI units.  The other code will
also perform similar unit conversions into its preferred units.

LAMMPS can also be used as an MDI driver in other unit choices it
supports, e.g. *lj*, but then no unit conversion is performed.

Related commands
""""""""""""""""

:doc:`mdi plugin <mdi>`, :doc:`mdi engine <mdi>`

Default
"""""""

The default for the optional keywords are virial = no, add = yes,
every = 1, connect = yes.
