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

At each contact, the elastic normal force of the *normal* model and the
normal damping force of the *damping* model act along the contact
normal at the contact point, scaled by the weight of the contact as for
:doc:`pair_style body/rounded/polyhedron <pair_body_rounded_polyhedron>`.
The cohesive force of :ref:`Wang <pair-gran-Wang>` with the strength
:math:`k_{na}` of the *cohesion* keyword is added to the normal force,
for surface separations up to the *cutoff*, and is scaled with the size
of the contact region, controlled by *delta_ua*.

As with the *history* keyword of :doc:`pair_style body/rounded/polyhedron
<pair_body_rounded_polyhedron>`, there is a single friction force per
pair of particles in contact, the tangential force of the *tangential*
model.  It acts at the average of the contact points, with the average
tangential relative velocity and contact normal, all weighted by the
elastic normal forces of the contacts.  The tangential displacement of
*linear_history* is stored per pair of particles.  The Coulomb limit of
the tangential force is :math:`\mu` times the sum of the elastic normal
forces of the contacts, without the damping forces, unlike :doc:`pair_style
granular <pair_granular>`, where the critical force includes the normal
damping.  The tangential damping is that of the *tangential* model, with
the damping prefactor of the *damping* model scaled by :math:`x_{\gamma,t}`.
The relative velocity is evaluated at the contact point between the
rounded surfaces, including the rotation of both particles.

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
