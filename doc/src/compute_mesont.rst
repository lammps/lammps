.. index:: compute mesont

compute mesont/Es command
==========================

compute mesont/Eb command
==========================

compute mesont/Et command
==========================

compute mesont/B  command
==========================

compute mesont/Es\_tot command
===============================

compute mesont/Eb\_tot command
===============================

compute mesont/Et\_tot command
===============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID mesont/Es

* ID, group-ID are documented in :doc:`compute <compute>` command
* mesont/Es = style name of the compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all mesont/Es

Description
"""""""""""

These computes define computations for the per-node stretching (mesont/Es),
bending (mesont/Eb), and intertube (mesont/Et) energies, buckling flag (mesont/B),
as well as the total stretching (mesont/Es\_tot), bending (mesont/Eb\_tot), and
intertube (mesont/Et\_tot) energies for each atom (node) in a group.

**Output info:**

These computes calculate per-node (per-atom) vectors (mesont/Es, mesont/Eb, mesont/Et, mesont/B), 
which can be accessed by any command that uses per-atom values from a 
compute as input, and global scalars (mesont/Es\_tot, mesont/Eb\_tot, 
mesont/Et\_tot). See the :doc:`Howto output <Howto_output>` doc page for an 
overview of LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
These computes are part of the USER-MESONT package. They are only enabled if 
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page for more info. In addition, :doc:`mesont pair_style <pair_style>`
must be used.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
