.. index:: fix smd/wall_surface

fix smd/wall_surface command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID smd/wall_surface arg type mol-ID

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/wall_surface = style name of this fix command
* arg = *file*

  .. parsed-literal::

        *file* = file name of a triangular mesh in stl format

* type = particle type to be given to the new particles created by this fix
* mol-ID = molecule-ID to be given to the new particles created by this fix (must be >= 65535)

Examples
""""""""

.. code-block:: LAMMPS

   fix stl_surf all smd/wall_surface tool.stl 2 65535

Description
"""""""""""

This fix creates reads a triangulated surface from a file in .STL
format.  For each triangle, a new particle is created which stores the
barycenter of the triangle and the vertex positions.  The radius of
the new particle is that of the minimum circle which encompasses the
triangle vertices.

The triangulated surface can be used as a complex rigid wall via the
:doc:`smd/tri_surface <pair_smd_triangulated_surface>` pair style.  It
is possible to move the triangulated surface via the
:doc:`smd/move_tri_surf <fix_smd_move_triangulated_surface>` fix style.

Immediately after a .STL file has been read, the simulation needs to
be run for 0 timesteps in order to properly register the new particles
in the system. See the "funnel_flow" example in the MACHDYN examples
directory.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to use Smooth Mach
Dynamics in LAMMPS.

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix *smd/wall\_surface* supports the *fix* keyword of :doc:`dump image
<dump_image>`.  The fix will pass geometry information about the wall
particles to *dump image* so that they be included in the rendered
image.

The color of the wall mesh object is by default that of the first atom
type when using color styles "type" or "element".  With color style
"const" the default value of "white" can be changed using
:doc:`dump_modify fcolor <dump_image>`.  The transparency is by default
fully opaque and can be changed with *dump\_modify ftrans*\ .

The *fflag1* setting of *dump image fix* determines whether the wall will
be rendered as a set of connected triangles (1) or as a mesh of cylinders (2).

When rendering triangles, the *fflag2* setting is ignored.  When using a
mesh of cylinders, the *fflag2* setting determines the diameter of the
cylinders.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Currently, no part of MACHDYN supports restarting nor
minimization. This fix has no outputs.

Restrictions
""""""""""""

This fix is part of the MACHDYN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The molecule ID given to the particles created by this fix have to be
equal to or larger than 65535.

Within each .STL file, only a single triangulated object must be
present, even though the STL format allows for the possibility of
multiple objects in one file.

Related commands
""""""""""""""""

:doc:`smd/triangle_mesh_vertices <compute_smd_triangle_vertices>`,
:doc:`smd/move_tri_surf <fix_smd_move_triangulated_surface>`,
:doc:`smd/tri_surface <pair_smd_triangulated_surface>`

Default
"""""""

none
