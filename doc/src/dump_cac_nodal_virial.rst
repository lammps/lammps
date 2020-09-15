.. index:: dump cac/nodal/virial

dump cac/nodal/virial command
=============================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/nodal/virial N file

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/nodal/virial = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to

Examples
""""""""

.. code-block:: LAMMPS

   dump mesh all cac/nodal/virial 1000 CACmesh.txt

Description
"""""""""""

Periodically outputs the equivalent nodal values of the virial stress. The nodal virial
of all finite elements for each of their respective internal degrees of freedom will be 
printed; atom positions will also be printed and denoted by the 0 element type following
one line of position data.

.. note::

   The element header for each atom or finite element possesses three additional
   dummy zeroes after the element scales due to implementation constraints. Make sure
   to read these in first in any post-processing script to avoid erroneous nodal info.
   The six vector of stress components is output in the following order: xx, yy, zz,
   xy, xz, yz

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.
Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/nodal/positions <dump_cac_initial_nodes>`,
:doc:`dump cac/kinetic/energy <dump_cac_kinetic>`, :doc:`dump cac/xyz <dump_cac_xyz>`,
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
