.. index:: atom_style

atom_style command
==================

Syntax
""""""


.. code-block:: LAMMPS

   atom_style style args

* style = *angle* or *atomic* or *body* or *bond* or *charge* or *dipole* or         *dpd* or *edpd* or *mdpd* or *tdpd* or *electron* or *ellipsoid* or         *full* or *line* or *meso* or *molecular* or *peri* or *smd* or         *sphere* or *spin* or *tri* or *template* or *hybrid*
  
  .. parsed-literal::
  
       args = none for any style except the following
         *body* args = bstyle bstyle-args
           bstyle = style of body particles
           bstyle-args = additional arguments specific to the bstyle
                         see the :doc:`Howto body <Howto_body>` doc page for details
         *tdpd* arg = Nspecies
           Nspecies = # of chemical species
         *template* arg = template-ID
           template-ID = ID of molecule template specified in a separate :doc:`molecule <molecule>` command
         *hybrid* args = list of one or more sub-styles, each with their args

* accelerated styles (with same args) = *angle/kk* or *atomic/kk* or *bond/kk* or *charge/kk* or *full/kk* or *molecular/kk*


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
   atom_style tdpd 2

Description
"""""""""""

Define what style of atoms to use in a simulation.  This determines
what attributes are associated with the atoms.  This command must be
used before a simulation is setup via a :doc:`read_data <read_data>`,
:doc:`read_restart <read_restart>`, or :doc:`create_box <create_box>`
command.

.. note::

   Many of the atom styles discussed here are only enabled if
   LAMMPS was built with a specific package, as listed below in the
   Restrictions section.

Once a style is assigned, it cannot be changed, so use a style general
enough to encompass all attributes.  E.g. with style *bond*\ , angular
terms cannot be used or added later to the model.  It is OK to use a
style more general than needed, though it may be slightly inefficient.

The choice of style affects what quantities are stored by each atom,
what quantities are communicated between processors to enable forces
to be computed, and what quantities are listed in the data file read
by the :doc:`read_data <read_data>` command.

These are the additional attributes of each style and the typical
kinds of physical systems they are used to model.  All styles store
coordinates, velocities, atom IDs and types.  See the
:doc:`read_data <read_data>`, :doc:`create_atoms <create_atoms>`, and
:doc:`set <set>` commands for info on how to set these various
quantities.

+--------------+-----------------------------------------------------+--------------------------------------+
| *angle*      | bonds and angles                                    | bead-spring polymers with stiffness  |
+--------------+-----------------------------------------------------+--------------------------------------+
| *atomic*     | only the default values                             | coarse-grain liquids, solids, metals |
+--------------+-----------------------------------------------------+--------------------------------------+
| *body*       | mass, inertia moments, quaternion, angular momentum | arbitrary bodies                     |
+--------------+-----------------------------------------------------+--------------------------------------+
| *bond*       | bonds                                               | bead-spring polymers                 |
+--------------+-----------------------------------------------------+--------------------------------------+
| *charge*     | charge                                              | atomic system with charges           |
+--------------+-----------------------------------------------------+--------------------------------------+
| *dipole*     | charge and dipole moment                            | system with dipolar particles        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *dpd*        | internal temperature and internal energies          | DPD particles                        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *edpd*       | temperature and heat capacity                       | eDPD particles                       |
+--------------+-----------------------------------------------------+--------------------------------------+
| *mdpd*       | density                                             | mDPD particles                       |
+--------------+-----------------------------------------------------+--------------------------------------+
| *tdpd*       | chemical concentration                              | tDPD particles                       |
+--------------+-----------------------------------------------------+--------------------------------------+
| *electron*   | charge and spin and eradius                         | electronic force field               |
+--------------+-----------------------------------------------------+--------------------------------------+
| *ellipsoid*  | shape, quaternion, angular momentum                 | aspherical particles                 |
+--------------+-----------------------------------------------------+--------------------------------------+
| *full*       | molecular + charge                                  | bio-molecules                        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *line*       | end points, angular velocity                        | rigid bodies                         |
+--------------+-----------------------------------------------------+--------------------------------------+
| *meso*       | rho, e, cv                                          | SPH particles                        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *mesont*     | mass, radius, length, buckling, connections, tube id| mesoscopic nanotubes                 |
+--------------+-----------------------------------------------------+--------------------------------------+
| *molecular*  | bonds, angles, dihedrals, impropers                 | uncharged molecules                  |
+--------------+-----------------------------------------------------+--------------------------------------+
| *peri*       | mass, volume                                        | mesoscopic Peridynamic models        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *smd*        | volume, kernel diameter, contact radius, mass       | solid and fluid SPH particles        |
+--------------+-----------------------------------------------------+--------------------------------------+
| *sphere*     | diameter, mass, angular velocity                    | granular models                      |
+--------------+-----------------------------------------------------+--------------------------------------+
| *spin*       | magnetic moment                                     | system with magnetic particles       |
+--------------+-----------------------------------------------------+--------------------------------------+
| *template*   | template index, template atom                       | small molecules with fixed topology  |
+--------------+-----------------------------------------------------+--------------------------------------+
| *tri*        | corner points, angular momentum                     | rigid bodies                         |
+--------------+-----------------------------------------------------+--------------------------------------+
| *wavepacket* | charge, spin, eradius, etag, cs\_re, cs\_im         | AWPMD                                |
+--------------+-----------------------------------------------------+--------------------------------------+

.. note::

   It is possible to add some attributes, such as a molecule ID, to
   atom styles that do not have them via the :doc:`fix property/atom <fix_property_atom>` command.  This command also
   allows new custom attributes consisting of extra integer or
   floating-point values to be added to atoms.  See the :doc:`fix property/atom <fix_property_atom>` doc page for examples of cases
   where this is useful and details on how to initialize, access, and
   output the custom values.

All of the above styles define point particles, except the *sphere*\ ,
*ellipsoid*\ , *electron*\ , *peri*\ , *wavepacket*\ , *line*\ , *tri*\ , and
*body* styles, which define finite-size particles.  See the :doc:`Howto spherical <Howto_spherical>` doc page for an overview of using
finite-size particle models with LAMMPS.

All of the point-particle styles assign mass to particles on a
per-type basis, using the :doc:`mass <mass>` command, The finite-size
particle styles assign mass to individual particles on a per-particle
basis.

For the *sphere* style, the particles are spheres and each stores a
per-particle diameter and mass.  If the diameter > 0.0, the particle
is a finite-size sphere.  If the diameter = 0.0, it is a point
particle.  Note that by use of the *disc* keyword with the :doc:`fix nve/sphere <fix_nve_sphere>`, :doc:`fix nvt/sphere <fix_nvt_sphere>`,
:doc:`fix nph/sphere <fix_nph_sphere>`, :doc:`fix npt/sphere <fix_npt_sphere>` commands, spheres can be effectively
treated as 2d discs for a 2d simulation if desired.  See also the :doc:`set density/disc <set>` command.

For the *ellipsoid* style, the particles are ellipsoids and each
stores a flag which indicates whether it is a finite-size ellipsoid or
a point particle.  If it is an ellipsoid, it also stores a shape
vector with the 3 diameters of the ellipsoid and a quaternion 4-vector
with its orientation.

For the *dipole* style, a point dipole is defined for each point
particle.  Note that if you wish the particles to be finite-size
spheres as in a Stockmayer potential for a dipolar fluid, so that the
particles can rotate due to dipole-dipole interactions, then you need
to use atom\_style hybrid sphere dipole, which will assign both a
diameter and dipole moment to each particle.

For the *electron* style, the particles representing electrons are 3d
Gaussians with a specified position and bandwidth or uncertainty in
position, which is represented by the eradius = electron size.

For the *peri* style, the particles are spherical and each stores a
per-particle mass and volume.

The *dpd* style is for dissipative particle dynamics (DPD) particles.
Note that it is part of the USER-DPD package, and is not for use with
the :doc:`pair_style dpd or dpd/stat <pair_dpd>` commands, which can
simply use atom\_style atomic.  Atom\_style dpd extends DPD particle
properties with internal temperature (dpdTheta), internal conductive
energy (uCond), internal mechanical energy (uMech), and internal
chemical energy (uChem).

The *edpd* style is for energy-conserving dissipative particle
dynamics (eDPD) particles which store a temperature (edpd\_temp), and
heat capacity(edpd\_cv).

The *mdpd* style is for many-body dissipative particle dynamics (mDPD)
particles which store a density (rho) for considering
density-dependent many-body interactions.

The *tdpd* style is for transport dissipative particle dynamics (tDPD)
particles which store a set of chemical concentration. An integer
"cc\_species" is required to specify the number of chemical species
involved in a tDPD system.

The *meso* style is for smoothed particle hydrodynamics (SPH)
particles which store a density (rho), energy (e), and heat capacity
(cv).

The *smd* style is for a general formulation of Smooth Particle
Hydrodynamics.  Both fluids and solids can be modeled.  Particles
store the mass and volume of an integration point, a kernel diameter
used for calculating the field variables (e.g. stress and deformation)
and a contact radius for calculating repulsive forces which prevent
individual physical bodies from penetrating each other.

For the *spin* style, a magnetic spin is associated to each atom.
Those spins have a norm (their magnetic moment) and a direction.

The *wavepacket* style is similar to *electron*\ , but the electrons may
consist of several Gaussian wave packets, summed up with coefficients
cs= (cs\_re,cs\_im).  Each of the wave packets is treated as a separate
particle in LAMMPS, wave packets belonging to the same electron must
have identical *etag* values.

For the *line* style, the particles are idealized line segments and
each stores a per-particle mass and length and orientation (i.e. the
end points of the line segment).

For the *tri* style, the particles are planar triangles and each
stores a per-particle mass and size and orientation (i.e. the corner
points of the triangle).

For the *mesont* style, the particles represent nodes of Nanotube
segments, and each stores a per-particle mass, radius, segment
length, tube id, buckling flag, and connections with neighbor nodes.

The *template* style allows molecular topology (bonds,angles,etc) to be
defined via a molecule template using the :doc:`molecule <molecule>`
command.  The template stores one or more molecules with a single copy
of the topology info (bonds,angles,etc) of each.  Individual atoms
only store a template index and template atom to identify which
molecule and which atom-within-the-molecule they represent.  Using the
*template* style instead of the *bond*\ , *angle*\ , *molecular* styles
can save memory for systems comprised of a large number of small
molecules, all of a single type (or small number of types).  See the
paper by Grime and Voth, in :ref:`(Grime) <Grime>`, for examples of how this
can be advantageous for large-scale coarse-grained systems.

.. note::

   When using the *template* style with a :doc:`molecule template <molecule>` that contains multiple molecules, you should
   insure the atom types, bond types, angle\_types, etc in all the
   molecules are consistent.  E.g. if one molecule represents H2O and
   another CO2, then you probably do not want each molecule file to
   define 2 atom types and a single bond type, because they will conflict
   with each other when a mixture system of H2O and CO2 molecules is
   defined, e.g. by the :doc:`read_data <read_data>` command.  Rather the
   H2O molecule should define atom types 1 and 2, and bond type 1.  And
   the CO2 molecule should define atom types 3 and 4 (or atom types 3 and
   2 if a single oxygen type is desired), and bond type 2.

For the *body* style, the particles are arbitrary bodies with internal
attributes defined by the "style" of the bodies, which is specified by
the *bstyle* argument.  Body particles can represent complex entities,
such as surface meshes of discrete points, collections of
sub-particles, deformable objects, etc.

The :doc:`Howto body <Howto_body>` doc page describes the body styles
LAMMPS currently supports, and provides more details as to the kind of
body particles they represent.  For all styles, each body particle
stores moments of inertia and a quaternion 4-vector, so that its
orientation and position can be time integrated due to forces and
torques.

Note that there may be additional arguments required along with the
*bstyle* specification, in the atom\_style body command.  These
arguments are described on the :doc:`Howto body <Howto_body>` doc page.


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
rotate due to torque, you need to use "atom\_style hybrid sphere
dipole".  When a hybrid style is used, atoms store and communicate the
union of all quantities implied by the individual styles.

When using the *hybrid* style, you cannot combine the *template* style
with another molecular style that stores bond,angle,etc info on a
per-atom basis.

LAMMPS can be extended with new atom styles as well as new body
styles; see the :doc:`Modify <Modify>` doc page.


----------


Styles with a *kk* suffix are functionally the same as the
corresponding style without the suffix.  They have been optimized to
run faster, depending on your available hardware, as discussed in on
the :doc:`Speed packages <Speed_packages>` doc page.  The accelerated
styles take the same arguments and should produce the same results,
except for round-off and precision issues.

Note that other acceleration packages in LAMMPS, specifically the GPU,
USER-INTEL, USER-OMP, and OPT packages do not use accelerated atom
styles.

The accelerated styles are part of the KOKKOS package.  They are only
enabled if LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

Restrictions
""""""""""""


This command cannot be used after the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command.

Many of the styles listed above are only enabled if LAMMPS was built
with a specific package, as listed below.  See the :doc:`Build package <Build_package>` doc page for more info.

The *angle*\ , *bond*\ , *full*\ , *molecular*\ , and *template* styles are
part of the MOLECULE package.

The *line* and *tri* styles are part of the ASPHERE package.

The *body* style is part of the BODY package.

The *dipole* style is part of the DIPOLE package.

The *peri* style is part of the PERI package for Peridynamics.

The *electron* style is part of the USER-EFF package for :doc:`electronic force fields <pair_eff>`.

The *dpd* style is part of the USER-DPD package for dissipative
particle dynamics (DPD).

The *edpd*\ , *mdpd*\ , and *tdpd* styles are part of the USER-MESO package
for energy-conserving dissipative particle dynamics (eDPD), many-body
dissipative particle dynamics (mDPD), and transport dissipative particle
dynamics (tDPD), respectively.

The *meso* style is part of the USER-SPH package for smoothed particle
hydrodynamics (SPH).  See `this PDF guide <USER/sph/SPH_LAMMPS_userguide.pdf>`_ to using SPH in LAMMPS.

The *mesont* style is part of the USER-MESONT package.

The *spin* style is part of the SPIN package.

The *wavepacket* style is part of the USER-AWPMD package for the
:doc:`antisymmetrized wave packet MD method <pair_awpmd>`.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`pair_style <pair_style>`

Default
"""""""

atom\_style atomic


----------


.. _Grime:



**(Grime)** Grime and Voth, to appear in J Chem Theory & Computation
(2014).
