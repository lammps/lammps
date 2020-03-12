.. index:: fix_modify AtC mesh write

fix_modify AtC mesh write command
===================================

Syntax
""""""

.. parsed-literal::

   fix_modify <AtC fixID> mesh write <f|p> <f|p> <f|p>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mesh write = name of the AtC sub-command
* filename = name of the file containing the mesh to be write

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mesh write myMesh.mesh

Description
"""""""""""

Writes a mesh to a text file.

Restrictions
""""""""""""

None

Related AtC commands
""""""""""""""""""""

:doc:`fix_modify AtC mesh create <atc_mesh_create>`
:doc:`fix_modify AtC mesh read <atc_mesh_read>`

Default
"""""""

None.
