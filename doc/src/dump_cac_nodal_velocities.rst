.. index:: dump cac/nodal/velocities

dump cac/nodal/velocities command
=================================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/nodal/velocities N file

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* cac/nodal/velocities = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to

Examples
""""""""

.. code-block:: LAMMPS

   dump meshvel all cac/nodal/velocities 1000 CACmeshvel.txt

Description
"""""""""""

Periodically outputs the nodal velocities at the specified interval. The nodal velocities
of all finite elements for each of their respective internal degrees of freedom will be 
printed; atom velocities will also be printed and denoted by the 0 element type following
one line of velocity data.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.

Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/initial/nodes <dump_cac_initial_nodes>`,
:doc:`dump cac/kinetic/energy <dump_cac_kinetic>`, :doc:`dump cac/xyz <dump_cac_xyz>`,
:doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
