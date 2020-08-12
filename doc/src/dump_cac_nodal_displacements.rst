.. index:: dump cac/nodal/displacements

dump cac/nodal/displacements command
====================================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/nodal/displacements N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/nodal/positions = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = same as arguments for :doc:`dump_style custom <dump>`

Examples
""""""""

.. code-block:: LAMMPS
   dump mesh all cac/nodal/displacements 1000 CACdis.txt

Description
"""""""""""

Periodically outputs the nodal displacements at the specified interval. The nodal displacements
of all finite elements for each of their respective internal degrees of freedom will be 
printed; atom displacements will also be printed and denoted by the 0 element type following
one line of displacement data. By default the displacement is computed with respect to the
reference coordinates established in the last read of a data file. To reset the reference coordinates
use the :doc:`reset_initial_nodes <reset_initial_nodes>` command.


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
