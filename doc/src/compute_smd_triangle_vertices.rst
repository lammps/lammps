.. index:: compute smd/triangle/vertices

compute smd/triangle/vertices command
=====================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/triangle/vertices

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/triangle/vertices = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/triangle/mesh/vertices

Description
"""""""""""

Define a computation that returns the coordinates of the vertices
corresponding to the triangle-elements of a mesh created by the :doc:`fix smd/wall_surface <fix_smd_wall_surface>`.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute returns a per-particle vector of vectors, which can be
accessed by any command that uses per-particle values from a compute
as input. See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

The per-particle vector has nine entries, (x1/y1/z1), (x2/y2/z2), and
(x3/y3/z3) corresponding to the first, second, and third vertex of
each triangle.

It is only meaningful to use this compute for a group of particles
which is created via the :doc:`fix smd/wall_surface <fix_smd_wall_surface>` command.

The output of this compute can be used with the dump2vtk_tris tool to
generate a VTK representation of the smd/wall_surface mesh for
visualization purposes.

The values will be given in :doc:`units <units>` of distance.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix smd/move/tri/surf <fix_smd_move_triangulated_surface>`,
:doc:`fix smd/wall_surface <fix_smd_wall_surface>`

Default
"""""""

none
