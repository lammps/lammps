.. index:: compute cnt

compute cnt/Es command
=======================

compute cnt/Eb command
=======================

compute cnt/Et command
=======================

compute cnt/B  command
=======================

compute cnt/Es\_tot command
============================

compute cnt/Eb\_tot command
============================

compute cnt/Et\_tot command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID cnt/Es

* ID, group-ID are documented in :doc:`compute <compute>` command
* cnt/Es = style name of the compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all cnt/Es

Description
"""""""""""

These computes define computations for the per-node stretching (cnt/Es),
bending (cnt/Eb), and intertube (cnt/Et) energies, buckling flag (cnt/B),
as well as the total stretching (cnt/Es\_tot), bending (cnt/Eb\_tot), and
intertube (cnt/Et\_tot) energies for each atom (node) in a group.

**Output info:**

These computes calculate per-node (per-atom) vectors (cnt/Es, cnt/Eb, cnt/Et, cnt/B), 
which can be accessed by any command that uses per-atom values from a 
compute as input, and global scalars (cnt/Es\_tot, cnt/Eb\_tot, 
cnt/Et\_tot). See the :doc:`Howto output <Howto_output>` doc page for an 
overview of LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
These computes are part of the USER-CNT package. They are only enabled if 
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page for more info. In addition, :doc:`cnt pair_style <pair_style>`
must be used.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
