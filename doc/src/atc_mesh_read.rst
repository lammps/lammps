.. index:: fix_modify AtC mesh read

fix_modify AtC mesh read command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> mesh read <f|p> <f|p> <f|p>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* mesh read = name of the AtC sub-command
* filename = name of the file containing the mesh to be read
* f or p = periodicity flags for x-, y-, and z-direction (optional)

Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC mesh read myComponent.mesh p p p
   fix_modify AtC mesh read myOtherComponent.exo

Description
"""""""""""

Reads a mesh from a text or exodus file, and assigns periodic boundary
conditions if needed.

Restrictions
""""""""""""

None

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`
- :doc:`fix_modify AtC mesh create <atc_mesh_create>`
- :doc:`fix_modify AtC mesh write <atc_mesh_write>`

Default
"""""""

Periodicity flags are set to false (f) by default.
