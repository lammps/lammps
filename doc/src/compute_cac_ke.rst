.. index:: compute cac/ke

compute cac/ke command
======================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp

* ID, group-ID are documented in :doc:`compute <compute>` command
* cac/ke = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all cac/ke
   compute myke mobile cac/ke

Description
"""""""""""

Define a computation that calculates the kinetic energy of a group of
atoms/elements by using their velocities and nodal velocities.
The kinetic energy of a finite element is estimated using the average
(with respect to node count) of the nodal velocities squared multiplied
by the element scales.

The number of atoms/elements contributing to the kinetic energy is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute\_modify <compute_modify>` command if this is not the case.

----------

**Output info:**

This compute calculates a global scalar (the kinetic energy).
This value can be used by any command that uses a global scalar 
as input. See the :doc:`Howto output <Howto_output>` doc page for
an overview of LAMMPS output options.

The scalar value calculated by this compute is "intensive".
The scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute requires a CAC atom style

**Default:** none
