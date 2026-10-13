.. index:: pair_style body/rounded/polyhedron/granular
.. index:: pair_style body/rounded/polyhedron/granular/omp

pair_style body/rounded/polyhedron/granular command
===================================================

Accelerator Variants: *body/rounded/polyhedron/granular/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style body/rounded/polyhedron/granular delta_ua cutoff

.. parsed-literal::

   delta_ua = multiple contact scaling factor
   cutoff = global separation cutoff for the cohesive interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style body/rounded/polyhedron/granular 1.0 0.5
   pair_coeff * * hooke 100.0 20.0 damping velocity tangential linear_history 30.0 0.25 0.5 cohesion 1.0

   pair_style body/rounded/polyhedron/granular 1.0 0.5
   pair_coeff 1 1 hertz/material 1200.0 20.0 0.3 damping mass_velocity &
              tangential linear_history 30.0 0.25 0.5 cohesion 1.0
   pair_coeff 2 2 hertz/material 800.0 20.0 0.25 damping mass_velocity &
              tangential linear_history 30.0 0.25 0.5 cohesion 1.0

Description
"""""""""""

.. versionadded:: TBD

Style *body/rounded/polyhedron/granular* computes the interactions
between 3d body particles of style *rounded/polyhedron* with the
normal, damping, and tangential contact models of :doc:`pair_style
granular <pair_granular>`.  The particles, their contacts, the cohesive
force, and the scaling of the cohesive force with the size of the contact
region are those of :doc:`pair_style body/rounded/polyhedron
<pair_body_rounded_polyhedron>` with the *history* keyword, whose page
describes them in detail.  Only the forces at the contacts are computed
differently, so that e.g. Hertzian contacts with the elastic properties
of the particles can be used instead of the linear springs of
:ref:`Wang <pair-gran-Wang>`.

The polyhedra are rounded at their vertices and edges by spheres with
the diameter specified in the data file.  Two particles interact at one
or more contacts, between a vertex, an edge, or a face of one particle
and a vertex, an edge, or a face of the other particle, with the overlap
:math:`\delta = -\delta_n` of their rounded surfaces at that contact.
For the normal contact models, each contact is treated as a contact
between two curved surfaces with the effective radius of curvature
:math:`R_{eff}`, which depends on the features in contact, where a
rounded vertex is a sphere and a rounded edge is a cylinder, both with
the rounded radius :math:`r` of the particle, and a face is flat:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Contact
     - :math:`R_{eff}`
   * - vertex - vertex
     - :math:`r_i r_j / (r_i + r_j)`
   * - vertex - face
     - :math:`r_v`, the radius of the particle with the vertex
   * - vertex - edge
     - :math:`r_v \sqrt{r_e / (r_v + r_e)}`
   * - edge - edge
     - :math:`\sqrt{r_i r_j}`

Here, :math:`1/R_{eff}` is the geometric mean of the principal relative
curvatures of the two surfaces, so that the elliptic contact region has
the same area as the circular one of a sphere with the radius
:math:`R_{eff}` on a flat surface :ref:`(Johnson) <pair-gran-Johnson>`.
For two edges, the curvatures are those of two cylinders crossing at a
right angle.  Two nearly parallel edges touch along a line, which is
represented by the contacts at the vertices at its ends instead, see
:doc:`pair_style body/rounded/polyhedron <pair_body_rounded_polyhedron>`.
A sphere, i.e. a body particle with a single vertex, touches another
particle at a single vertex, edge, or face contact.

Since two polyhedra can touch at several points at once, the contact
models are applied at two levels: the normal force and the normal
damping act at each contact point, while the friction force acts once
per pair of particles, from all its contacts together.  For a pair of
particles, the forces are computed in these steps:

#. **Contacts.**  The contacts between the two particles are found as
   for :doc:`pair_style body/rounded/polyhedron
   <pair_body_rounded_polyhedron>`.  Each contact has a contact point
   on either particle, which defines the contact normal, the surface
   separation :math:`\delta_n`, the kind of features in contact, which
   defines :math:`R_{eff}` as in the table above, and a weight
   :math:`w \le 1`.  The weights are smaller than one for the corners of
   the region of overlap of nearly parallel faces and for the crossings
   of nearly parallel edges, so that these contacts appear and disappear
   smoothly as the particles move.  Contacts at the same point are
   counted once.

#. **Normal force at each contact.**  Each contact with an overlap
   :math:`\delta = -\delta_n > 0` is a separate contact of the *normal*
   model, with its own :math:`\delta` and :math:`R_{eff}`, and the
   effective mass of the two particles.  The force

   .. math::

      F_n = w \left( F_{ne}(\delta, R_{eff}) + j_a F_{c} \right)

   acts along the contact normal at the contact point, where
   :math:`F_{ne}` is the elastic force of the *normal* model and
   :math:`F_c = -k_{na} (r_c - \delta_n)` is the cohesive force of
   :ref:`Wang <pair-gran-Wang>` with the strength :math:`k_{na}` of the
   *cohesion* keyword and the *cutoff* :math:`r_c`.  The cohesive force
   is scaled by the factor :math:`j_a \ge 1`, which grows with the size of
   the contact region, controlled by *delta_ua*, and is reduced near the
   ends of crossing edges, as described for :doc:`pair_style
   body/rounded/polyhedron <pair_body_rounded_polyhedron>`.  For surface
   separations :math:`0 < \delta_n \le r_c`, only the cohesive force acts.
   The total normal force of two particles touching at several points is
   the sum of these forces.  E.g. a cube resting with a face on another
   cube touches it at the corners of the region of overlap, each of them
   with its own Hertzian force, rather than as one flat contact, so that
   the stiffness of such a contact depends on the number and the weights
   of these corners.

#. **Normal damping at each contact.**  At each contact in overlap, the
   *damping* model adds the force :math:`w F_{damp}` along the contact
   normal, from the normal relative velocity of the two particles at the
   contact point, including their rotation.  There is no tangential
   damping at the individual contacts.

#. **Friction force of the pair of particles.**  As with the *history*
   keyword of :doc:`pair_style body/rounded/polyhedron
   <pair_body_rounded_polyhedron>`, there is a single friction force per
   pair of particles in contact.  The contacts in overlap are averaged,
   each weighted by its elastic normal force :math:`w F_{ne}`: the contact
   points give the point where the friction force acts, and the contact
   normals and tangential relative velocities give the normal and the
   tangential velocity for the *tangential* model.  A single friction
   force avoids jumps of the friction force between contacts, e.g. when
   the contact with the largest overlap changes.  The tangential
   displacement of *linear_history* is stored per pair of particles and
   is reset when they no longer touch.  The Coulomb limit of the
   tangential force is :math:`\mu \sum w F_{ne}`, the friction
   coefficient times the sum of the elastic normal forces of the contacts,
   without the damping forces, unlike :doc:`pair_style granular
   <pair_granular>`, where the critical force includes the normal damping.
   The tangential damping is that of the *tangential* model, with the
   damping prefactor of the *damping* model scaled by
   :math:`x_{\gamma,t}`.  The friction force acts on both particles with
   opposite signs at the average contact point, so that it also exerts
   torques.

The contact point lies halfway between the two rounded surfaces in
overlap, on the contact normal, i.e. at the distance
:math:`r - \delta/2` from the vertex, from the nearest point on the edge,
or from the face, and not from the center of mass of the particle.  The
forces at a contact act at this point, so that their torque on a particle
is the cross product of the vector from its center of mass to the contact
point with the force.  Since the center of mass of a polyhedron is in
general not on the contact normal, also the normal forces exert torques.
The relative velocities are evaluated at the contact point as well,
including the rotation of both particles, while :doc:`pair_style granular
<pair_granular>` uses the distance :math:`r` from the center of a sphere
for the rotational part of the relative velocity.  For two spherical
body particles, i.e. with a single vertex, without rotation, the forces
are the same as those of :doc:`pair_style granular <pair_granular>` with
the same coefficients.

The coefficients of each pair of atom types are specified with the
:doc:`pair_coeff <pair_coeff>` command, with the same syntax for the
contact models as for :doc:`pair_style granular <pair_granular>`:

.. code-block:: LAMMPS

   pair_coeff I J normal_model normal_coeffs tangential tangential_model tangential_coeffs keyword values ...

* normal_model = *hooke* or *hertz* or *hertz/material*
* normal_coeffs = coefficients of the normal model, see :doc:`pair_style granular <pair_granular>`
* tangential_model = *linear_nohistory* or *linear_history*
* tangential_coeffs = coefficients of the tangential model, see :doc:`pair_style granular <pair_granular>`
* zero or more keyword/value pairs may be appended
* keyword = *damping* or *cohesion*

  .. parsed-literal::

       *damping* value = *velocity* or *mass_velocity*
       *cohesion* value = k_na
         k_na = strength of the cohesive force (energy/distance\^2 units)

The *tangential* keyword with its model is required, while the *damping*
and *cohesion* keywords are optional.  The meaning of the coefficients
is the same as for :doc:`pair_style granular <pair_granular>`, with the overlap and the effective radius of
curvature of each contact as described above.  Other normal, damping,
and tangential models, as well as the *rolling* and *twisting* models,
*limit_damping*, and *synchronized_verlet* of :doc:`pair_style granular
<pair_granular>` are not supported, since they are not well defined for
the multiple contacts of polyhedra or for a single friction force per pair
of particles.

The input script *in.heap3d_granular* in the *examples/body* directory
pours cubes and tetrahedra onto a floor with Hertzian contacts from the
elastic properties of either particle type, and tilts gravity once they
have come to rest in a heap.  The floor uses the same contact models with
:doc:`fix wall/body/polyhedron/granular <fix_wall_body_polyhedron_granular>`.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The coefficients of two different atom types I,J that are not specified
by a pair_coeff command are mixed from those of I,I and J,J as for
:doc:`pair_style granular <pair_granular>`, and the strength of the
cohesive force as the geometric mean.  This pair style does not support
the :doc:`pair_modify <pair_modify>` shift, table, and tail options.

The total energy is not conserved, as described for :doc:`pair_style
body/rounded/polyhedron <pair_body_rounded_polyhedron>`.  The reported
pair energy includes the elastic energy of the normal models and the
energy of the cohesive force.  As for that pair style, the work done by
the forces that do not derive from this energy can be accessed by the
:doc:`compute pair <compute_pair>` command.

This pair style does not write its information to :doc:`binary restart
files <restart>`.  Thus, you need to re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.  The
tangential displacements of the contacts are written to the restart file
and are used again when the pair style is re-specified.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

This pair style is part of the BODY package and requires the GRANULAR
package.  It is only enabled if LAMMPS was built with both packages.  See
the :doc:`Build package <Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style body/rounded/polyhedron
<pair_body_rounded_polyhedron>`, :doc:`pair_style granular
<pair_granular>`, :doc:`fix wall/body/polyhedron
<fix_wall_body_polyhedron>`, :doc:`fix wall/body/polyhedron/granular
<fix_wall_body_polyhedron_granular>`

Default
"""""""

The *damping* model is *velocity*, and the strength of the cohesive force
is :math:`k_{na} = 0`.

----------

.. _pair-gran-Wang:

**(Wang)** J. Wang, H. S. Yu, P. A. Langston, F. Y. Fraige, Granular
Matter, 13, 1 (2011).

.. _pair-gran-Johnson:

**(Johnson)** K. L. Johnson, Contact Mechanics, Cambridge University
Press (1985).
