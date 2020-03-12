.. index:: compute pe

compute pe command
==================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID pe keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* pe = style name of this compute command
* zero or more keywords may be appended
* keyword = *pair* or *bond* or *angle* or *dihedral* or *improper* or *kspace* or *fix*

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pe
   compute molPE all pe bond angle dihedral improper

Description
"""""""""""

Define a computation that calculates the potential energy of the
entire system of atoms.  The specified group must be "all".  See the
:doc:`compute pe/atom <compute_pe_atom>` command if you want per-atom
energies.  These per-atom values could be summed for a group of atoms
via the :doc:`compute reduce <compute_reduce>` command.

The energy is calculated by the various pair, bond, etc potentials
defined for the simulation.  If no extra keywords are listed, then the
potential energy is the sum of pair, bond, angle, dihedral, improper,
kspace (long-range), and fix energy.  I.e. it is as if all the
keywords were listed.  If any extra keywords are listed, then only
those components are summed to compute the potential energy.

The Kspace contribution requires 1 extra FFT each timestep the energy
is calculated, if using the PPPM solver via the :doc:`kspace_style pppm <kspace_style>` command.  Thus it can increase the cost of the
PPPM calculation if it is needed on a large fraction of the simulation
timesteps.

Various fixes can contribute to the total potential energy of the
system if the *fix* contribution is included.  See the doc pages for
:doc:`individual fixes <fix>` for details of which ones compute a
potential energy.

.. note::

   The :doc:`fix_modify energy yes <fix_modify>` command must also be
   specified if a fix is to contribute potential energy to this command.

A compute of this style with the ID of "thermo\_pe" is created when
LAMMPS starts up, as if this command were in the input script:

.. code-block:: LAMMPS

   compute thermo_pe all pe

See the "thermo\_style" command for more details.

----------

**Output info:**

This compute calculates a global scalar (the potential energy).  This
value can be used by any command that uses a global scalar value from
a compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute pe/atom <compute_pe_atom>`

**Default:** none
