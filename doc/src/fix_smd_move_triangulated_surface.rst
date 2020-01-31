.. index:: fix smd/move\_tri\_surf

fix smd/move\_tri\_surf command
===============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID smd/move_tri_surf keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/move\_tri\_surf keyword = style name of this fix command
* keyword = *\*LINEAR* or *\*WIGGLE* or *\*ROTATE*
  
  .. parsed-literal::
  
        *\*LINEAR* args = Vx Vy Vz
           Vx,Vy,Vz = components of velocity vector (velocity units), any component can be specified as NULL
        *\*WIGGLE* args = Vx Vy Vz max_travel
           vx,vy,vz = components of velocity vector (velocity units), any component can be specified as NULL
           max_travel = wiggle amplitude
        *\*ROTATE* args = Px Py Pz Rx Ry Rz period
           Px,Py,Pz = origin point of axis of rotation (distance units)
           Rx,Ry,Rz = axis of rotation vector
           period = period of rotation (time units)



Examples
""""""""


.. parsed-literal::

   fix 1 tool smd/move_tri_surf \*LINEAR 20 20 10
   fix 2 tool smd/move_tri_surf \*WIGGLE 20 20 10
   fix 2 tool smd/move_tri_surf \*ROTATE 0 0 0 5 2 1

Description
"""""""""""

This fix applies only to rigid surfaces read from .STL files via fix
:doc:`smd/wall\_surface <fix_smd_wall_surface>` .  It updates position
and velocity for the particles in the group each timestep without
regard to forces on the particles.  The rigid surfaces can thus be
moved along simple trajectories during the simulation.

The *\*LINEAR* style moves particles with the specified constant velocity
vector V = (Vx,Vy,Vz). This style also sets the velocity of each particle
to V = (Vx,Vy,Vz).

The *\*WIGGLE* style moves particles in an oscillatory fashion.
Particles are moved along (vx, vy, vz) with constant velocity until a
displacement of max\_travel is reached. Then, the velocity vector is
reversed. This process is repeated.

The *\*ROTATE* style rotates particles around a rotation axis R =
(Rx,Ry,Rz) that goes through a point P = (Px,Py,Pz). The period of the
rotation is also specified. This style also sets the velocity of each
particle to (omega cross Rperp) where omega is its angular velocity
around the rotation axis and Rperp is a perpendicular vector from the
rotation axis to the particle.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth Mach
Dynamics in LAMMPS.

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no part of USER-SMD supports restarting nor
minimization. This fix has no outputs.

Restrictions
""""""""""""


This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`smd/triangle\_mesh\_vertices <compute_smd_triangle_vertices>`,
:doc:`smd/wall\_surface <fix_smd_wall_surface>`

**Default:** none


