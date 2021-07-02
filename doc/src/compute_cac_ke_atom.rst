.. index:: compute cac/ke/atom

compute cac/ke/atom command
===========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp

* ID, group-ID are documented in :doc:`compute <compute>` command
* cac/ke/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all cac/ke/atom
   compute myke mobile cac/ke/atom

Description
"""""""""""

Define a computation that calculates the kinetic energy of each atom or
element in the specified group of atoms/elements by using their velocities
and nodal velocities. The kinetic energy of a finite element is estimated 
using the average (with respect to node count) of the nodal kinetic energies 
multiplied by the element scales.

The value of the kinetic energy will be 0.0 for atoms and elements
not in the specified compute group.

----------

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute requires a CAC atom style

**Default:** none
