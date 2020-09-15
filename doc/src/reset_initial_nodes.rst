.. index:: reset_initial_nodes

reset_initial_nodes command
=============================

Syntax
""""""

.. parsed-literal::

   reset_initial_nodes

Examples
""""""""

.. code-block:: LAMMPS

   reset_initial_nodes

Description
"""""""""""

Reset the initial nodal positions to the current nodal nodal positions. This is useful
in the context of computing nodal displacements with :doc:`dump cac/nodal/displacements <dump_cac_nodal_displacements>`
so that the displacement is computed with respect to an updated set of reference coordinates. This
command requires a cac atom style and should be used after reading a data or restart file.

Restrictions
""""""""""""
 none

**Default:** none
