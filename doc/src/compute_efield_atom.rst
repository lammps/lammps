.. index:: compute efield/atom

compute efield/atom command
===========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID efield/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* efield/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all efield/atom
   compute 1 all efield/atom pair yes kspace no

Used in input scripts:

   .. parsed-literal::

      examples/PACKAGES/dielectric/in.confined
      examples/PACKAGES/dielectric/in.nopbc

Description
"""""""""""

Define a computation that calculates the electric field at each atom in a group.
The compute should only enabled with pair and kspace styles that are provided
by the DIELECTRIC package because only these styles compute the per-atom
electric field at every time step.

The electric field is a 3-component vector.  The value of the electric field
components will be 0.0 for atoms not in the specified compute group.

----------

The keyword/value option pairs are used in the following ways.

For the *pair* and *kspace* keywords, the real-space and reciprocal-space
contributions to the electric field can be turned off and on.


Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values will be in electric field :doc:`units <units>`.

Restrictions
""""""""""""
This compute is part of the DIELECTRIC package. It is only enabled if
LAMMPS was built with that package.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

Default
"""""""

The option defaults are pair = yes and kspace = yes.

