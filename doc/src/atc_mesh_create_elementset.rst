.. index:: fix_modify AtC mesh create_elementset

fix_modify AtC mesh create_elementset command
=============================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> mesh create_elementset <id> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mesh create_elementset = name of the AtC sub-command
* id = id to assign to the collection of FE nodes
* <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of the bounding box that contains only the desired elements

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mesh create_elementset middle -4.1 4.1 -100 100 -100 1100

Description
"""""""""""

Command to assign an id to a set of FE elements to be used subsequently
in defining material and mesh-based operations.

Restrictions
""""""""""""

Only viable for rectangular grids.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC mesh delete_elements <atc_mesh_delete_elements>`
- :doc:`fix_modify AtC mesh nodeset_to_elementset <atc_mesh_nodeset_to_elementset>`

Default
"""""""

Coordinates are assumed to be in lattice units.
