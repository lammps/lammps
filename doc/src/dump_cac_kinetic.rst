.. index:: dump cac/kinetic/energy

dump cac/kinetic/energy command
===============================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/kinetic/energy N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/kinetic/energy = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = same as arguments for :doc:`dump_style custom <dump>`

Examples
""""""""

.. code-block:: LAMMPS

   dump meshvel all cac/kinetic/energy 1000 CACmeshvel.txt

Description
"""""""""""

Periodically outputs the per node kinetic energy due to each dimension of velocity
at the specified interval. This is printed for all finite elements for each of their respective
internal degrees of freedom; for atoms, denoted by the 0 element type, there will be one line
of data after the one line declaring the atom's type information.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.

Restrictions
""""""""""""

This dump style requires a CAC :doc:`atom style <atom_style>`

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/initial/nodes <dump_cac_initial_nodes>`,
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`, :doc:`dump cac/xyz <dump_cac_xyz>`,
:doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
