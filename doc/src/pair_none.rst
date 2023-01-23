.. index:: pair_style none

pair_style none command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style none

Examples
""""""""

.. code-block:: LAMMPS

   pair_style none

Description
"""""""""""

Using a pair style of *none* means that any previous pair style setting
will be deleted and pairwise forces and energies are not computed.

As a consequence there will be a pairwise force cutoff of 0.0, which has
implications for the default setting of the neighbor list and the
communication cutoff.  Those are the sum of the largest pairwise cutoff
and the neighbor skin distance (see the documentation of the
:doc:`neighbor <neighbor>` command and the :doc:`comm_modify
<comm_modify>` command).  When you have bonds, angles, dihedrals, or
impropers defined at the same time, you must set the communication
cutoff so that communication cutoff distance is large enough to acquire
and communicate sufficient ghost atoms from neighboring subdomains as
needed for computing bonds, angles, etc.

A pair style of *none* will also not request a pairwise neighbor list.
However if the :doc:`neighbor <neighbor>` style is *bin*, data
structures for binning are still allocated.  If the neighbor list cutoff
is small, then these data structures can consume a large amount of
memory.  So you should either set the neighbor style to *nsq* or set the
skin distance to a larger value.

See the :doc:`pair_style zero <pair_zero>` for a way to set a pairwise
cutoff and thus trigger the building of a neighbor lists and setting
a corresponding communication cutoff, but compute no pairwise interactions.

Restrictions
""""""""""""

You must not use a :doc:`pair_coeff <pair_coeff>` command with this pair
style.  Since there is no interaction computed, you cannot set any
coefficients for it.

Related commands
""""""""""""""""

:doc:`pair_style zero <pair_zero>`

Default
"""""""

none
