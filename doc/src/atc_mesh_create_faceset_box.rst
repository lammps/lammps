.. index:: fix_modify AtC mesh create_faceset box

fix_modify AtC mesh create_faceset box command
==============================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> mesh create_faceset <id> box <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <in|out> [units]

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mesh create_faceset = name of the AtC sub-command
* id = id to assign to the collection of FE faces
* box = use bounding box to define FE faces
* <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of the bounding box that is coincident with the desired FE faces
* <in|out> = "in" gives inner faces to the box, "out" gives the outer faces to the box
* units = option to specify real as opposed to lattice units


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mesh create_faceset obndy box -4.0 4.0 -12 12 -12 12 out


Description
"""""""""""

Command to assign an id to a set of FE faces.

Restrictions
""""""""""""

Only viable for rectangular grids.

Related AtC commands
""""""""""""""""""""

:doc:`fix_modify AtC mesh create_faceset plane <atc_mesh_create_faceset_plane>`

Default
"""""""

The default options are units = lattice and the use of outer faces.
