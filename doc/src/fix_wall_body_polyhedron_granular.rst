.. index:: fix wall/body/polyhedron/granular

fix wall/body/polyhedron/granular command
=========================================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/body/polyhedron/granular normal_model normal_coeffs tangential tangential_model tangential_coeffs keyword values ... wallstyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/body/polyhedron/granular = style name of this fix command
* normal_model = *hooke* or *hertz* or *hertz/material*
* normal_coeffs = coefficients of the normal model, see :doc:`pair_style granular <pair_granular>`
* tangential_model = *linear_nohistory* or *linear_history*
* tangential_coeffs = coefficients of the tangential model, see :doc:`pair_style granular <pair_granular>`
* zero or more keyword/value pairs may be appended to the models
* keyword = *damping*

  .. parsed-literal::

       *damping* value = *velocity* or *mass_velocity*

* wallstyle = *xplane* or *yplane* or *zplane*
* args = list of arguments for a particular style

  .. parsed-literal::

       *xplane* or *yplane* or *zplane* args = lo hi
         lo,hi = position of lower and upper plane (distance units), either can be NULL)

* zero or more keyword/value pairs may be appended to args
* keyword = *wiggle*

  .. parsed-literal::

       *wiggle* values = dim amplitude period
         dim = *x* or *y* or *z*
         amplitude = size of oscillation (distance units)
         period = time of oscillation (time units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all wall/body/polyhedron/granular hooke 2000.0 50.0 tangential linear_history 500.0 0.25 0.5 zplane 0.0 NULL
   fix 1 all wall/body/polyhedron/granular hertz/material 1200.0 20.0 0.3 damping mass_velocity &
         tangential linear_history 30.0 0.25 0.5 zplane 0.0 NULL wiggle x 0.5 10.0

Description
"""""""""""

.. versionadded:: TBD

This fix bounds the simulation domain of 3d body particles of style
*rounded/polyhedron* with planar walls, as :doc:`fix wall/body/polyhedron
<fix_wall_body_polyhedron>` does, but with the normal, damping, and
tangential contact models of :doc:`pair_style granular <pair_granular>`.
It is the counterpart of :doc:`pair_style body/rounded/polyhedron/granular
<pair_body_rounded_polyhedron_granular>` for the walls, e.g. to use
Hertzian contacts with the elastic properties of the particles for both
the contacts between particles and with the walls.

As for :doc:`fix wall/body/polyhedron <fix_wall_body_polyhedron>`, each
vertex of a particle, or the center of a sphere, touches the wall when
its distance from the wall is smaller than its rounded radius :math:`r`.
The contact is that of a sphere with the radius :math:`r` with a flat
wall, i.e. the effective radius of curvature is :math:`R_{eff} = r`, with
the overlap :math:`\delta = r - d`, where :math:`d` is the signed distance
of the vertex from the wall.  At each such contact, the elastic normal
force of the *normal* model and the normal damping force of the *damping*
model act along the wall normal at the contact point between the rounded
surface and the wall.  The wall has the same elastic properties as the
particles, as for :doc:`fix wall/gran <fix_wall_gran>`.  There is no
cohesion with the wall.

As with the *history* keyword of :doc:`fix wall/body/polyhedron
<fix_wall_body_polyhedron>`, the friction force of the wall on a particle
is a single force, the tangential force of the *tangential* model, which
acts at the average of the contact points of the vertices touching the
wall, with their average tangential velocity relative to the wall, both
weighted by their elastic normal forces.  Its Coulomb limit is
:math:`\mu` times the sum of the elastic normal forces of these vertices.
The tangential displacement of *linear_history* is stored per particle
and wall, and is reset to zero when the particle no longer touches the
wall.  The tangential damping is that of the *tangential* model, as for
:doc:`pair_style body/rounded/polyhedron/granular
<pair_body_rounded_polyhedron_granular>`.  The relative velocity and the
torques are evaluated at the contact point between the rounded surface and
the wall, at the distance :math:`r - \delta/2` from the vertex, while
:doc:`fix wall/gran <fix_wall_gran>` uses the radius of the particle for
both.

The *tangential* keyword with its model is required and the *damping*
model is optional.  Other normal, damping, and tangential models, as well
as the *rolling*, *twisting*, and *heat* models, *limit_damping*, and
*synchronized_verlet* of :doc:`pair_style granular <pair_granular>` are
not supported.  The *wallstyle* and the *wiggle* keyword are the same as
for :doc:`fix wall/body/polyhedron <fix_wall_body_polyhedron>`.

The input script *in.heap3d_granular* in the *examples/body* directory
uses this fix for the floor on which cubes and tetrahedra come to rest
in a heap.

-----------------

Dump image info
"""""""""""""""

This fix supports the *fix* keyword of :doc:`dump image <dump_image>`
in the same way as :doc:`fix wall/body/polyhedron
<fix_wall_body_polyhedron>`.

------------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the tangential displacements of the particles touching
the walls to :doc:`binary restart files <restart>`, so that a simulation
can continue correctly.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy
minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BODY package and requires the GRANULAR package.
It is only enabled if LAMMPS was built with both packages.  See the
:doc:`Build package <Build_package>` page for more info.

Any dimension (xyz) that has a wall must be non-periodic.

Related commands
""""""""""""""""

:doc:`fix wall/body/polyhedron <fix_wall_body_polyhedron>`,
:doc:`pair_style body/rounded/polyhedron/granular
<pair_body_rounded_polyhedron_granular>`, :doc:`fix wall/gran
<fix_wall_gran>`, :doc:`pair_style granular <pair_granular>`

Default
"""""""

The *damping* model is *velocity*.
