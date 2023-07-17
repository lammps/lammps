.. index:: compute rattlers

compute rattlers command
========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID rattlers cutoff zmin ntries

* ID, group-ID are documented in :doc:`compute <compute>` command
* rattlers = style name of this compute command
* cutoff = *type* or *radius*

  .. parsed-literal::

       *type* = cutoffs determined based on atom types
       *radius* = cutoffs determined based on atom diameters (atom style sphere)

* zmin = minimum coordination for a non-rattler particle
* ntries = maximum number of iterations to remove rattlers

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all rattlers type 4 10

Description
"""""""""""

Define a compute that identifies rattlers in a system. Rattlers are
identified using an interative approach. The coordination number of
all atoms is first calculated.  The *type* and *radius* settings are
used to select whether interaction cutoffs are determined by atom
types or by the sum of atomic radii (atom style sphere), respectively.
Rattlers are then identified as particles with a coordination number
less than *zmin* and are removed from consideration. Atomic coordination
numbers are then recalculated, excluding previously identified rattlers,
to identify a new set of rattlers. This process is iterated, up to a maximum
of *ntries*, until no new rattlers are identified and the remaining
atoms form a stable network of contacts.

Output info
"""""""""""

This compute calculates a per-atom vector and a global scalar. The vector
designates which atoms are rattlers, indicated by a value 1. Non-rattlers
have a value of 0. The global scalar returns the total number of rattlers
in the system. See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

The *radius* cutoff option requires that atoms store a radius as defined by the
:doc:`atom_style sphere <atom_style>` or similar commands.

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`
:doc:`compute contact/atom <compute_contact_atom>`

Default
"""""""

none
