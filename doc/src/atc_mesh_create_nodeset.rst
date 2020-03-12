.. index:: fix_modify AtC mesh create_nodeset

fix_modify AtC mesh create_nodeset command
==========================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> mesh create_nodeset <id> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mesh create\_nodeset = name of the AtC sub-command
* id = id to assign to the collection of FE nodes
* <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of the bounding box that contains only the desired nodes

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mesh create_nodeset lbc -12.1 -11.9 -12 12 -12 12

Description
"""""""""""

Command to assign an id to a set of FE nodes to be used subsequently in
defining boundary conditions.

Restrictions
""""""""""""

None

Related AtC commands
""""""""""""""""""""

:doc:`fix_modify AtC mesh add_to_nodeset <atc_mesh_add_to_nodeset>`

Default
"""""""

Coordinates are assumed to be in lattice units.
