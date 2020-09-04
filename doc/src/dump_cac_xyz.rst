.. index:: dump cac/xyz

dump cac/xyz command
====================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/xyz N file

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/xyz = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to

Examples
""""""""

.. code-block:: LAMMPS

   dump mesh all cac/xyz 1000 CACmesh.txt

Description
"""""""""""

Periodically outputs the position of a set of virtual atoms and real atoms for the 
specified group. These virtual atoms are the result of converting a CAC element into 
atomistic resolution through the interpolation fields. The purpose of this is to 
print small regions of a CAC model that may require the use of an atomistic 
diagnostic tool such as OVITO.

.. warning::

   CAC elements can be quite large, as a result a conversion to atomistic
   resolution may result in a massive amount of data if your group size is not carefully
   chosen. Currently you cannot elect to output a piece of a CAC element as virtual atoms,
   if an element is included in the group (because its centroid is in the group region etc.)
   then all of its material will be converted to virtual atoms and output.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.

Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/initial/nodes <dump_cac_initial_nodes>`,
:doc:`dump cac/atom <dump_cac_atom>`, :doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
