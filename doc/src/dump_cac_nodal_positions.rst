.. index:: dump cac/nodal/positions

dump cac/nodal/positions command
================================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/nodal/positions N file

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/nodal/positions = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to

Examples
""""""""

.. code-block:: LAMMPS

   dump mesh all cac/nodal/positions 1000 CACmesh.txt

Description
"""""""""""

Periodically outputs the nodal positions at the specified interval. The nodal positions
of all finite elements for each of their respective internal degrees of freedom will be 
printed; atom positions will also be printed and denoted by the 0 element type following
one line of position data.

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
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
