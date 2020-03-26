.. index:: compute mesont

compute mesont command
==========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID mesont mode

* ID, group-ID are documented in :doc:`compute <compute>` command
* mesont = style name of the compute command
* mode = one of estretch, ebend, etube, stretch_tot, ebend_tot, and etube_tot (see details below)

Examples
""""""""


.. parsed-literal::

   compute 1 all mesont estretch

Description
"""""""""""

These computes define computations for the per-node stretching (estretch),
bending (ebend), and intertube (etube) energies, as well as the total 
stretching (estretch_tot), bending (ebend_tot), and intertube (etube_tot) 
energies for each atom (node) in a group. The evaluated value is selected by 
a parameter passed to the compute: estretch, ebend, etube, estretch_tot, 
ebend_tot, and etube_tot.

**Output info:**

These computes calculate per-node (per-atom) vectors (estretch, ebend, etube), 
which can be accessed by any command that uses per-atom values from a 
compute as input, and global scalars (stretch_tot, ebend_tot, and etube_tot). 
See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS 
output options.

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
