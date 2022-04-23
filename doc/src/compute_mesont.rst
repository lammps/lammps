.. index:: compute mesont

compute mesont command
======================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID mesont mode

* ID, group-ID are documented in :doc:`compute <compute>` command
* mesont = style name of the compute command
* mode = one of estretch, ebend, etube (see details below)

Examples
""""""""


.. parsed-literal::

   compute 1 all mesont estretch

Description
"""""""""""

These computes define computations for the stretching (estretch), bending
(ebend), and intertube (etube) per-node (atom) and total energies. The
evaluated value is selected by a parameter passed to the compute: estretch,
ebend, etube.

Output info
"""""""""""

These computes calculate per-node (per-atom) vectors, which can be accessed by
any command that uses per-atom values from a compute as input, and global
scalars. See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The computed values are provided in energy :doc:`units <units>`.

Restrictions
""""""""""""
These computes are part of the MESONT package. They are only enabled if
LAMMPS is built with that package. See the :doc:`Build package <Build_package>`
doc page for more info. In addition, :doc:`mesont pair_style <pair_style>`
must be used.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

Default
"""""""

none

