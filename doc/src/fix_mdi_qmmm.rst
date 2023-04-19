.. index:: fix mdi/qmmm

fix mdi/qmmm command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mdi/qmmm mode keyword value(s) keyword value(s) ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mdi/qmmm = style name of this fix command
* mode = *direct* or *potential*
* zero or more keyword/value pairs may be appended
* keyword = *virial* or *add* or *every* or *connect* or *elements*

  .. parsed-literal::

       *virial* args = *yes* or *no*
         yes = request virial tensor from server code
         no = do not request virial tensor from server code
       *connect* args = *yes* or *no*
         yes = perform a one-time connection to the MDI engine code
         no = do not perform the connection operation
       *elements* args = N_1 N_2 ... N_ntypes
         N_1,N_2,...N_ntypes = chemical symbol for each of ntypes LAMMPS atom types

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mdi/qmmm direct
   fix 1 all mdi/qmmm potential virial yes
   fix 1 all mdi/qmmm potential virial yes elements 13 29

Description
"""""""""""

.. versionadded:: 28Mar2023

This command enables LAMMPS to act as a client with another server code
to perform a coupled QM/MM (quantum-mechanics/molecular-mechanics)
simulation.  LAMMPS will perform classical MD (molecular mechanics
or MM) for the (typically larger) MM portion of the system.  A quantum
mechanics code will calculate quantum energy and forces for the QM
portion of the system.  The two codes work together to calculate the
energy and forces due to the cross interactions between QM and MM atoms.
The QM server code must support use of the `MDI Library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_ as
explained below.

The partitioning of the system between QM and MM atoms is as follows.
Atoms in the specified group are QM atoms; the remaining atoms are MM
atoms.  The input script should thus define this partitioning.
See additional information below about other requirements for an input
script to use this fix and perform a QM/MM simulation.

The code coupling performed by this command is done via the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_.
LAMMPS runs as an MDI driver (client), and sends MDI commands to an
external MDI engine code (server), in this case a QM code which has
support for MDI.  See the :doc:`Howto mdi <Howto_mdi>` page for more
information about how LAMMPS can operate as either an MDI driver or
engine.

The ``examples/QUANTUM`` directory has sub-directories with example
input scripts using this fix in tandem with different QM codes.  The
README files in the sub-directories explain how to download and build
the various QM codes.  They also explain how to launch LAMMPS and the QM
code so that they communicate via the MDI library using either MPI or
sockets.  Any QM code that supports MDI could be used in addition to
those discussed in the sub-directories.  See the :doc:`Howto mdi
<Howto_mdi>` page for a current list (March 2022) of such QM codes.

Note that an engine code can support MDI in either or both of two modes.
It can be used as a stand-alone code, launched at the same time as
LAMMPS.  Or it can be used as a plugin library, which LAMMPS loads.  See
the :doc:`mdi plugin <mdi>` command for how to trigger LAMMPS to load a
plugin library.  The ``examples/QUANTUM`` sub-directory README files
explains how to launch the two codes in either mode.

----------

The *mode* setting determines which QM/MM coupling algorithm is used.
LAMMPS currently supports *direct* and *potential* algorithms, based
on the *mode* setting.  Both algorithms should give reasonably
accurate results, but some QM codes support only one of the two modes.
E.g. in the ``examples/QUANTUM`` directory, PySCF supports only *direct*,
NWChem supports only *potential*, and LATTE currently supports
neither, so it cannot be used for QM/MM simulations using this fix.

The *direct* option passes the coordinates and charges of each MM atom
to the quantum code, in addition to the coordinates of each QM atom.
The quantum code returns forces on each QM atom as well as forces on
each MM atom.  The latter is effectively the force on MM atoms due to
the QM atoms.

The input script for performing a *direct* mode QM/MM simulation should
do the following:

* delete all bonds (angles, dihedrals, etc) between QM atoms
* set the charge on each QM atom to zero
* define no bonds (angles, dihedrals, etc) which involve both QM and MM atoms
* define a force field (pair, bonds, angles, optional kspace) for the entire system

The first two bullet can be performed using the :doc:`delete_bonds
<delete_bonds>` and :doc:`set <set>` commands.

The third bullet is required to have a consistent model, but is not
checked by LAMMPS.

The fourth bullet implies that non-bonded non-Coulombic interactions
(e.g. van der Waals) between QM/QM and QM/MM pairs of atoms are
computed by LAMMPS.

See the ``examples/QUANTUM/PySCF/in.*`` files for examples of input
scripts for QM/MM simulations using the *direct* mode.

The *potential* option passes the coordinates of each QM atom and a
Coulomb potential for each QM atom to the quantum code.  The latter is
calculated by performing a Coulombics-only calculation for the entire
system, subtracting all QM/QM pairwise Coulombic terms, and dividing
the Coulomb energy on each QM atom by the charge of the QM atom.  The
potential value represents the Coulombic influence of all the MM atoms
on each QM atom.

The quantum code returns forces and charge on each QM atom.  The new
charges on the QM atom are used to re-calculate the MM force field,
resulting in altered forces on the MM atoms.

The input script for performing a *potential* mode QM/MM simulation
should do the following:

* delete all bonds (angles, dihedrals, etc) between QM atoms
* define a hybrid pair style which includes a Coulomb-only pair sub-style
* define no bonds (angles, dihedrals, etc) which involve both QM and MM atoms
* define a force field (pair, bonds, angles, optional kspace) for the entire system

The first operation can be performed using the :doc:`delete_bonds
<delete_bonds>` command.  See the ``examples/QUANTUM/NWChem/in.*`` files
for examples of how to do this.

The second operation is necessary so that this fix can calculate the
Coulomb potential for the QM atoms.

The third bullet is required to have a consistent model, but is not
checked by LAMMPS.

The fourth bullet implies that non-bonded non-Coulombic interactions
(e.g. van der Waals) between QM/QM and QM/MM pairs of atoms are computed
by LAMMPS.  However, some QM codes do not want the MM code (LAMMPS) to
compute QM/QM van der Waals interactions.  NWChem is an example.  In
this case, the coefficients for those interactions need to be turned
off, which typically requires the atom types for the QM atoms be
different than those for the MM atoms.

See the ``examples/QUANTUM/NWChem/in.*`` files for examples of input
scripts for QM/MM simulations using the *potential* mode.  Those scripts
also illustrate how to turn off QM/QM van der Waals interactions.

----------

The *virial* keyword setting of yes or no determines whether LAMMPS
will request the QM code to also compute and return the QM
contribution to a stress tensor for the system which LAMMPS will
convert to a 6-element symmetric virial tensor.

The *connect* keyword determines whether this fix performs a one-time
connection to the QM code.  The default is *yes*.  The only time a
*no* is needed is if this command is used multiple times in an input
script.  E.g. if it used inside a loop which also uses the :doc:`clear
<clear>` command to destroy the system (including this fix).  As
example would be a script which loop over a series of independent QM/MM
simulations, e.g. each with their own data file.  In this use case
*connect no* could be used along with the :doc:`mdi connect and exit
<mdi>` command to one-time initiate/terminate the connection outside
the loop.

The *elements* keyword allows specification of what element each
LAMMPS atom type corresponds to.  This is specified by the chemical
symbol of the element, e.g. C or Al or Si.  A symbol must be specified
for each of the ntypes LAMMPS atom types.  Multiple LAMMPS types can
represent the same element.  Ntypes is typically specified via the
:doc:`create_box <create_box>` command or in the data file read by the
:doc:`read_data <read_data>` command.

If this keyword is specified, then this fix will send the MDI
">ELEMENTS" command to the engine, to insure the two codes are
consistent in their definition of atomic species.  If this keyword is
not specified, then this fix will send the MDI >TYPES command to the
engine.  This is fine if both the LAMMPS driver and the MDI engine are
initialized so that the atom type values are consistent in both codes.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the potential energy computed by the QM code to the
global potential energy of the system as part of :doc:`thermodynamic
output <thermo_style>`.  The default setting for this fix is
:doc:`fix_modify energy yes <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by
this fix to add the contribution computed by the QM code to the global
pressure of the system as part of :doc:`thermodynamic output
<thermo_style>`.  The default setting for this fix is :doc:`fix_modify
virial yes <fix_modify>`.

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
accessed by various :doc:`output commands <Howto_output>`.  Note that
for *direct* mode this will be quantum forces on both QM and MM atoms.
For *potential* mode it will only be quantum forces on QM atoms; the
forces for MM atoms will be zero.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces computed by the QM code are used during an energy
minimization, invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the potential energy associated with the QM forces to
   be included in the total potential energy of the system (the
   quantity being minimized), you MUST not disable the
   :doc:`fix_modify <fix_modify>` *energy* option for this fix.


Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To use LAMMPS as an MDI driver in conjunction with other MDI-enabled
codes (MD or QM codes), the :doc:`units <units>` command should be
used to specify *real* or *metal* units.  This will ensure the correct
unit conversions between LAMMPS and MDI units.  The other code will
also perform similar unit conversions into its preferred units.

If this fix is used in conjuction with a QM code that does not support
periodic boundary conditions (more specifically, a QM code that does
not support the ``>CELL`` MDI command), the LAMMPS system must be
fully non-periodic.  I.e. no dimension of the system can be periodic.

Related commands
""""""""""""""""

:doc:`mdi plugin <mdi>`,
:doc:`mdi engine <mdi>`,
:doc:`fix mdi/qm <fix_mdi_qm>`

Default
"""""""

The default for the optional keywords are virial = no and connect = yes.
