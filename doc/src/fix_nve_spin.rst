.. index:: fix nve/spin

fix nve/spin command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nve/spin keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/spin = style name of this fix command
* keyword = *lattice*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate both spin and atomic degress of freedom
         frozen = integrate spins on a fixed lattice

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all nve/spin lattice moving
   fix 1 all nve/spin lattice frozen

Description
"""""""""""

Perform a symplectic integration for the spin or spin-lattice system.

The *lattice* keyword defines if the spins are integrated on a lattice
of fixed atoms (lattice = frozen), or if atoms are moving
(lattice = moving).
The first case corresponds to a spin dynamics calculation, and
the second to a spin-lattice calculation.
By default a spin-lattice integration is performed (lattice = moving).

The *nve/spin* fix applies a Suzuki-Trotter decomposition to
the equations of motion of the spin lattice system, following the scheme:

.. image:: JPG/fix_integration_spin_stdecomposition.jpg
   :align: center

according to the implementation reported in :ref:`(Omelyan) <Omelyan1>`.

A sectoring method enables this scheme for parallel calculations.
The implementation of this sectoring algorithm is reported
in :ref:`(Tranchida) <Tranchida1>`.

----------

Restrictions
""""""""""""

This fix style can only be used if LAMMPS was built with the SPIN
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

To use the spin algorithm, it is necessary to define a map with
the atom\_modify command. Typically, by adding the command:

.. code-block:: LAMMPS

   atom_modify map array

before you create the simulation box. Note that the keyword "hash"
instead of "array" is also valid.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`fix nve <fix_nve>`

Default
"""""""

The option default is lattice = moving.

----------

.. _Omelyan1:

**(Omelyan)** Omelyan, Mryglod, and Folk. Phys. Rev. Lett.
86(5), 898. (2001).

.. _Tranchida1:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
