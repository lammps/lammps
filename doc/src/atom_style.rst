.. index:: atom_style

atom_style command
==================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style style args

* style = *amoeba* or *angle* or *atomic* or *body* or *bond* or *charge* or *dielectric* or *dipole* or  *dpd* or *edpd* or *electron* or *ellipsoid* or *full* or *line* or *mdpd* or *molecular* or *oxdna* or *peri* or *smd* or *sph* or *sphere* or *bpm/sphere* or *spin* or *tdpd* or *tri* or *template* or *wavepacket* or *hybrid*

  .. parsed-literal::

       args = none for any style except the following
         *body* args = bstyle bstyle-args
           bstyle = style of body particles
           bstyle-args = additional arguments specific to the bstyle
                         see the :doc:`Howto body <Howto_body>` doc
                         page for details
         *sphere* arg = 0/1 (optional) for static/dynamic particle radii
         *bpm/sphere* arg = 0/1 (optional) for static/dynamic particle radii
         *tdpd* arg = Nspecies
           Nspecies = # of chemical species
         *template* arg = template-ID
           template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
         *hybrid* args = list of one or more sub-styles, each with their args

* accelerated styles (with same args) = *angle/kk* or *atomic/kk* or *bond/kk* or *charge/kk* or *full/kk* or *molecular/kk* or *spin/kk*

Examples
""""""""

.. code-block:: LAMMPS

   atom_style atomic
   atom_style bond
   atom_style full
   atom_style body nparticle 2 10
   atom_style hybrid charge bond
   atom_style hybrid charge body nparticle 2 5
   atom_style spin
   atom_style template myMols
   atom_style hybrid template twomols charge
   atom_style tdpd 2

Description
"""""""""""

The *atom_style* command selects which per-atom attributes are
associated with atoms in a LAMMPS simulation and thus stored and
communicated with those atoms as well as read from and stored in data
and restart files.  Different models (e.g. :doc:`pair styles
<pair_style>`) require access to specific per-atom attributes and thus
require a specific atom style.  For example, to compute Coulomb
interactions, the atom must have a "charge" (aka "q") attribute.

A number of distinct atom styles exist that combine attributes.  Some
atom styles are a superset of other atom styles.  Further attributes
may be added to atoms either via using a hybrid style which provides a
union of the attributes of the sub-styles, or via the :doc:`fix
property/atom <fix_property_atom>` command.  The *atom_style* command
must be used before a simulation is setup via a :doc:`read_data
<read_data>`, :doc:`read_restart <read_restart>`, or :doc:`create_box
<create_box>` command.

.. note::

   Many of the atom styles discussed here are only enabled if LAMMPS was
   built with a specific package, as listed below in the Restrictions
   section.

Once a style is selected and the simulation box defined, it cannot be
changed but only augmented with the :doc:`fix property/atom
<fix_property_atom>` command.  So one should select an atom style
general enough to encompass all attributes required.  E.g. with atom
style *bond*, it is not possible to define angles and use angle styles.

It is OK to use a style more general than needed, though it may be
slightly inefficient because it will allocate and communicate
additional unused data.

Atom style attributes
"""""""""""""""""""""

The atom style *atomic* has the minimum subset of per-atom attributes
and is also the default setting.  It encompasses the following per-atom
attributes (name of the vector or array in the :doc:`Atom class
<Classes_atom>` is given in parenthesis): atom-ID (tag), type (type),
position (x), velocities (v), forces (f), image flags (image), group
membership (mask).  Since all atom styles are a superset of atom style
*atomic*\ , they all include these attributes.

This table lists all the available atom styles, which attributes they
provide, which :doc:`package <Packages>` is required to use them, and
what the typical applications are that use them.  See the
:doc:`read_data <read_data>`, :doc:`create_atoms <create_atoms>`, and
:doc:`set <set>` commands for details on how to set these various
quantities.  More information about many of the styles is provided in
the Additional Information section below.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Atom style
     - Attributes
     - Required package
     - Applications
   * - *amoeba*
     - *full* + "1-5 special neighbor data"
     - :ref:`AMOEBA <PKG-AMOEBA>`
     - AMOEBA/HIPPO force fields
   * - *angle*
     - *bond* + "angle data"
     - :ref:`MOLECULE <PKG-MOLECULE>`
     - bead-spring polymers with stiffness
   * - *atomic*
     - tag, type, x, v, f, image, mask
     -
     - atomic liquids, solids, metals
   * - *body*
     - *atomic* + radius, rmass, angmom, torque, body
     - :ref:`BODY <PKG-BODY>`
     - arbitrary bodies, see :doc:`body howto <Howto_body>`
   * - *bond*
     - *atomic* + molecule, nspecial, special + "bond data"
     - :ref:`MOLECULE <PKG-MOLECULE>`
     - bead-spring polymers
   * - *bpm/sphere*
     - *bond* + radius, rmass, omega, torque, quat
     - :ref:`BPM <PKG-BPM>`
     - granular bonded particle models, see :doc:`BPM howto <Howto_bpm>`
   * - *charge*
     - *atomic* + q
     -
     - atomic systems with charges
   * - *dielectric*
     - *full* + mu, area, ed, em, epsilon, curvature, q_scaled
     - :ref:`DIELECTRIC <PKG-DIELECTRIC>`
     - systems with surface polarization
   * - *dipole*
     - *charge* + mu
     - :ref:`DIPOLE <PKG-DIPOLE>`
     - atomic systems with charges and point dipoles
   * - *dpd*
     - *atomic* + rho + "reactive DPD data"
     - :ref:`DPD-REACT <PKG-DPD-REACT>`
     - reactive DPD
   * - *edpd*
     - *atomic* + "eDPD data"
     - :ref:`DPD-MESO <PKG-DPD-MESO>`
     - Energy conservative DPD (eDPD)
   * - *electron*
     - *charge* + espin, eradius, ervel, erforce
     - :ref:`EFF <PKG-EFF>`
     - Electron force field systems
   * - *ellipsoid*
     - *atomic* + rmass, angmom, torque, ellipsoid
     -
     - aspherical particles
   * - *full*
     - *molecular* + q
     - :ref:`MOLECULE <PKG-MOLECULE>`
     - molecular force fields
   * - *line*
     - *atomic* + molecule, radius, rmass, omega, torque, line
     -
     - 2-d rigid body particles
   * - *mdpd*
     - *atomic* + rho, drho, vest
     - :ref:`DPD-MESO <PKG-DPD-MESO>`
     - Many-body DPD (mDPD)
   * - *molecular*
     - *angle* + "dihedral and improper data"
     - :ref:`MOLECULE <PKG-MOLECULE>`
     - apolar and uncharged molecules
   * - *oxdna*
     - *atomic* + id5p
     - :ref:`CG-DNA <PKG-CG-DNA>`
     - coarse-grained DNA and RNA models
   * - *peri*
     - *atomic* + rmass, vfrac, s0, x0
     - :ref:`PERI <PKG-PERI>`
     - mesoscopic Peridynamics models
   * - *smd*
     - *atomic* + molecule, radius, rmass + "smd data"
     - :ref:`MACHDYN <PKG-MACHDYN>`
     - Smooth Mach Dynamics models
   * - *sph*
     - *atomic* + "sph data"
     - :ref:`SPH <PKG-SPH>`
     - Smoothed particle hydrodynamics models
   * - *sphere*
     - *atomic* + radius, rmass, omega, torque
     -
     - finite size spherical particles, e.g. granular models
   * - *spin*
     - *atomic* + "magnetic moment data"
     - :ref:`SPIN <PKG-SPIN>`
     - magnetic particles
   * - *tdpd*
     - *atomic* + cc, cc_flux, vest
     - :ref:`DPD-MESO <PKG-DPD-MESO>`
     - Transport DPD (tDPD)
   * - *template*
     - *atomic* + molecule, molindex, molatom
     - :ref:`MOLECULE <PKG-MOLECULE>`
     - molecular systems where attributes are taken from :doc:`molecule files <molecule>`
   * - *tri*
     - *sphere* + molecule, angmom, tri
     -
     - 3-d triangulated rigid body LJ particles
   * - *wavepacket*
     - *charge* + "wavepacket data"
     - :ref:`AWPMD <PKG-AWPMD>`
     - Antisymmetrized wave packet MD

.. note::

   It is possible to add some attributes, such as a molecule ID and
   charge, to atom styles that do not have them built in using the
   :doc:`fix property/atom <fix_property_atom>` command.  This command
   also allows new custom-named attributes consisting of extra integer
   or floating-point values or vectors to be added to atoms.  See the
   :doc:`fix property/atom <fix_property_atom>` page for examples of
   cases where this is useful and details on how to initialize,
   access, and output these custom values.

----------

Particle size and mass
""""""""""""""""""""""

All of the atom styles define point particles unless they (1) define
finite-size spherical particles via the *radius* attribute, or (2)
define finite-size aspherical particles (e.g. the *body*, *ellipsoid*,
*line*, and *tri* styles).  Most of these styles can also be used with
mixtures of point and finite-size particles.

Note that the *radius* property may need to be provided as a
*diameter* (e.g. in :doc:`molecule files <molecule>` or :doc:`data
files <read_data>`).  See the :doc:`Howto spherical <Howto_spherical>`
page for an overview of using finite-size spherical and aspherical
particle models with LAMMPS.

Unless an atom style defines the per-atom *rmass* attribute, particle
masses are defined on a per-type basis, using the :doc:`mass <mass>`
command.  This means each particle's mass is indexed by its atom
*type*.

A few styles define the per-atom *rmass* attribute which can also be
added using the :doc:`fix property/atom <fix_property_atom>` command.
In this case each particle stores its own mass.  Atom styles that have
a per-atom rmass may define it indirectly through setting particle
diameter and density on a per-particle basis.  If both per-type mass
and per-atom *rmass* are defined (e.g. in a hybrid style), the
per-atom mass will take precedence in any operation which which works
with both flavors of mass.

----------

Additional information about specific atom styles
"""""""""""""""""""""""""""""""""""""""""""""""""

For the *body* style, the particles are arbitrary bodies with internal
attributes defined by the "style" of the bodies, which is specified by
the *bstyle* argument.  Body particles can represent complex entities,
such as surface meshes of discrete points, collections of
sub-particles, deformable objects, etc.

The :doc:`Howto body <Howto_body>` page describes the body styles
LAMMPS currently supports, and provides more details as to the kind of
body particles they represent.  For all styles, each body particle
stores moments of inertia and a quaternion 4-vector, so that its
orientation and position can be time integrated due to forces and
torques.

Note that there may be additional arguments required along with the
*bstyle* specification, in the atom_style body command.  These
arguments are described on the :doc:`Howto body <Howto_body>` doc page.

For the *dielectric* style, each particle can be either a physical
particle (e.g. an ion), or an interface particle representing a boundary
element between two regions of different dielectric constant. For
interface particles, in addition to the properties associated with
atom_style full, each particle also should be assigned a normal unit
vector (defined by normx, normy, normz), an area (area/patch), the
difference and mean of the dielectric constants of two sides of the
interface along the direction of the normal vector (ed and em), the
local dielectric constant at the boundary element (epsilon), and a mean
local curvature (curv).  Physical particles must be assigned these
values, as well, but only their local dielectric constants will be used;
see documentation for associated :doc:`pair styles <pair_dielectric>`
and :doc:`fixes <fix_polarize>`.  The distinction between the physical
and interface particles is only meaningful when :doc:`fix polarize
<fix_polarize>` commands are applied to the interface particles. This
style is part of the DIELECTRIC package.

For the *dipole* style, a point dipole vector mu is defined for each
point particle.  Note that if you wish the particles to be finite-size
spheres as in a Stockmayer potential for a dipolar fluid, so that the
particles can rotate due to dipole-dipole interactions, then you need
to use the command `atom_style hybrid sphere dipole`, which will
assign both a diameter and dipole moment to each particle.  This also
requires using an integrator with a "/sphere" suffix like :doc:`fix
nve/sphere <fix_nve_sphere>` or :doc:`fix nvt/sphere <fix_nvt_sphere>`
and the "update dipole" or "update dlm" parameters to the fix
commands.

The *dpd* style is for reactive dissipative particle dynamics (DPD)
particles.  Note that it is part of the DPD-REACT package, and is not
required for use with the :doc:`pair_style dpd or dpd/stat <pair_dpd>`
commands, which only require the attributes from atom_style *atomic*.
Atom_style *dpd* extends DPD particle properties with internal
temperature (dpdTheta), internal conductive energy (uCond), internal
mechanical energy (uMech), and internal chemical energy (uChem).

The *edpd* style is for energy-conserving dissipative particle
dynamics (eDPD) particles which store a temperature (edpd_temp), and
heat capacity (edpd_cv).

For the *electron* style, the particles representing electrons are 3d
Gaussians with a specified position and bandwidth or uncertainty in
position, which is represented by the eradius = electron size.

For the *ellipsoid* style, particles can be ellipsoids which each
stores a shape vector with the 3 diameters of the ellipsoid and a
quaternion 4-vector with its orientation.  Each particle stores a flag
in the ellipsoid vector which indicates whether it is an ellipsoid (1)
or a point particle (0).

For the *line* style, particles can be are idealized line segments
which store a per-particle mass and length and orientation (i.e. the
end points of the line segment).  Each particle stores a flag in the
line vector which indicates whether it is a line segment (1) or a
point particle (0).

The *mdpd* style is for many-body dissipative particle dynamics (mDPD)
particles which store a density (rho) for considering density-dependent
many-body interactions.

The *oxdna* style is for coarse-grained nucleotides and stores the
3'-to-5' polarity of the nucleotide strand, which is set through
the bond topology in the data file. The first (second) atom in a
bond definition is understood to point towards the 3'-end (5'-end)
of the strand.

For the *peri* style, the particles are spherical and each stores a
per-particle mass and volume.

The *smd* style is for Smooth Particle Mach dynamics.  Both fluids and
solids can be modeled.  Particles store the mass and volume of an
integration point, a kernel diameter used for calculating the field
variables (e.g. stress and deformation) and a contact radius for
calculating repulsive forces which prevent individual physical bodies
from penetrating each other.

The *sph* style is for smoothed particle hydrodynamics (SPH) particles
which store a density (rho), energy (esph), and heat capacity (cv).

For the *spin* style, a magnetic spin is associated with each atom.
Those spins have a norm (their magnetic moment) and a direction.

The *tdpd* style is for transport dissipative particle dynamics (tDPD)
particles which store a set of chemical concentration. An integer
"cc_species" is required to specify the number of chemical species
involved in a tDPD system.

The *wavepacket* style is similar to the *electron* style, but the
electrons may consist of several Gaussian wave packets, summed up with
coefficients cs= (cs_re,cs_im).  Each of the wave packets is treated
as a separate particle in LAMMPS, wave packets belonging to the same
electron must have identical *etag* values.

The *sphere* and *bpm/sphere* styles allow particles to be either point
particles or finite-size particles.  If the *radius* attribute is >
0.0, the particle is a finite-size sphere.  If the diameter = 0.0, it
is a point particle.  Note that by using the *disc* keyword with the
:doc:`fix nve/sphere <fix_nve_sphere>`, :doc:`fix nvt/sphere
<fix_nvt_sphere>`, :doc:`fix nph/sphere <fix_nph_sphere>`, :doc:`fix
npt/sphere <fix_npt_sphere>` commands for the *sphere* style, spheres
can be effectively treated as 2d discs for a 2d simulation if desired.
See also the :doc:`set density/disc <set>` command.  These styles also
take an optional 0 or 1 argument.  A value of 0 means the radius of
each sphere is constant for the duration of the simulation (this is
the default).  A value of 1 means the radii may vary dynamically
during the simulation, e.g. due to use of the :doc:`fix adapt
<fix_adapt>` command.

The *template* style allows molecular topology (bonds,angles,etc) to be
defined via a molecule template using the :doc:`molecule <molecule>`
command.  The template stores one or more molecules with a single copy
of the topology info (bonds,angles,etc) of each.  Individual atoms only
store a template index and template atom to identify which molecule and
which atom-within-the-molecule they represent.  Using the *template*
style instead of the *bond*, *angle*, *molecular* styles can save memory
for systems comprised of a large number of small molecules, all of a
single type (or small number of types).  See the paper by Grime and
Voth, in :ref:`(Grime) <Grime>`, for examples of how this can be
advantageous for large-scale coarse-grained systems.  The
``examples/template`` directory has a few demo inputs and examples
showing the use of the *template* atom style versus *molecular*.

.. note::

   When using the *template* style with a :doc:`molecule template
   <molecule>` that contains multiple molecules, you should ensure the
   atom types, bond types, angle_types, etc in all the molecules are
   consistent.  E.g. if one molecule represents H2O and another CO2,
   then you probably do not want each molecule file to define 2 atom
   types and a single bond type, because they will conflict with each
   other when a mixture system of H2O and CO2 molecules is defined,
   e.g. by the :doc:`read_data <read_data>` command.  Rather the H2O
   molecule should define atom types 1 and 2, and bond type 1.  And
   the CO2 molecule should define atom types 3 and 4 (or atom types 3
   and 2 if a single oxygen type is desired), and bond type 2.

For the *tri* style, particles can be planar triangles which each
stores a per-particle mass and size and orientation (i.e. the corner
points of the triangle).  Each particle stores a flag in the tri
vector which indicates whether it is a triangle (1) or a point
particle (0).

----------

Typically, simulations require only a single (non-hybrid) atom style.
If some atoms in the simulation do not have all the properties defined
by a particular style, use the simplest style that defines all the
needed properties by any atom.  For example, if some atoms in a
simulation are charged, but others are not, use the *charge* style.
If some atoms have bonds, but others do not, use the *bond* style.

The only scenario where the *hybrid* style is needed is if there is no
single style which defines all needed properties of all atoms.  For
example, as mentioned above, if you want dipolar particles which will
rotate due to torque, you need to use "atom_style hybrid sphere
dipole".  When a hybrid style is used, atoms store and communicate the
union of all quantities implied by the individual styles.

When using the *hybrid* style, you cannot combine the *template* style
with another molecular style that stores bond, angle, etc info on a
per-atom basis.

LAMMPS can be extended with new atom styles as well as new body styles;
see the corresponding manual page on :doc:`modifying & extending LAMMPS
<Modify_atom>`.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command.

Many of the styles listed above are only enabled if LAMMPS was built
with a specific package, as listed below.  See the :doc:`Build package
<Build_package>` page for more info.  The table above lists which package
is required for individual atom styles.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`pair_style <pair_style>`,
:doc:`fix property/atom <fix_property_atom>`, :doc:`set <set>`

Default
"""""""

The default atom style is *atomic*.  If atom_style *sphere* or
*bpm/sphere* is used, its default argument is 0.

----------

.. _Grime:

**(Grime)** Grime and Voth, to appear in J Chem Theory & Computation
(2014).
