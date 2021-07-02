.. index:: dump cac/initial/nodes

dump cac/initial/nodes command
==============================

Syntax
""""""

.. parsed-literal::

   dump ID group-ID cac/initial/nodes N file args

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms/elements to be dumped
* cac/initial/nodes = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = same as arguments for :doc:`dump_style custom <dump>`

Examples
""""""""

.. code-block:: LAMMPS

   dump initmesh all cac/initial/nodes 1 initial_CACmesh.txt

Description
"""""""""""

Generates a dump file containing the initial nodal positions of the simulated system. This is 
typically required as reference information to plot deformation contours and more. Since the information
is unchanging this is simply a mechanism to output the information in case the initial data file is 
not available or impractical; as such it only needs to be dumped once.

.. warning::

   When restarting a CAC simulation the initial nodal
   positions will not be nodal positions at the beginning of the restart run. 
   They will be the nodal positions from the initial data file that was read or created
   for that model in the first run.

----------

NOTE: The :doc:`dump_modify sort <dump_modify>` option
does not work with this dump style.
Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump cac/nodal/positions <dump_cac_nodal_positions>`,
:doc:`dump cac/kinetic/energy <dump_cac_kinetic>`, :doc:`dump cac/xyz <dump_cac_xyz>`,
:doc:`dump cac/nodal/velocities <dump_cac_nodal_velocities>`,
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`
