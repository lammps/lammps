.. index:: set

set command
===========

Syntax
""""""

.. code-block:: LAMMPS

   set style ID keyword values ...

* style = *atom* or *type* or *mol* or *group* or *region*
* ID = depends on style

.. parsed-literal::

       for style = *atom*, ID = a range of atom IDs
       for style = *type*, ID = a range of numeric types or a single type label
       for style = *mol*, ID = a range of molecule IDs
       for style = *group*, ID = a group ID
       for style = *region*, ID = a region ID

* one or more keyword/value pairs may be appended
* keyword = *type* or *type/fraction* or *type/ratio* or *type/subset*
  or *mol* or *x* or *y* or *z* or *vx* or *vy* or *vz* or *charge* or
  *dipole* or *dipole/random* or *quat* or *spin/atom* or *spin/atom/random* or
  *spin/electron* or *radius/electron* or
  *quat* or *quat/random* or *diameter* or *shape* or *length* or *tri* or
  *theta* or *theta/random* or *angmom* or *omega* or
  *mass* or *density* or *density/disc* or *temperature* or
  *volume* or *image* or *bond* or *angle* or *dihedral* or
  *improper* or *sph/e* or *sph/cv* or *sph/rho* or
  *smd/contact/radius* or *smd/mass/density* or *dpd/theta* or
  *edpd/temp* or *edpd/cv* or *cc* or *epsilon* or
  *i_name* or *d_name* or *i2_name* or *d2_name*

  .. parsed-literal::

       *type* value = numeric atom type or type label
         value can be an atom-style variable (see below)
       *type/fraction* values = type fraction seed
         type = numeric atom type or type label
         fraction = approximate fraction of selected atoms to set to new atom type
         seed = random # seed (positive integer)
       *type/ratio* values = type fraction seed
         type = numeric atom type or type label
         fraction = exact fraction of selected atoms to set to new atom type
         seed = random # seed (positive integer)
       *type/subset* values = type Nsubset seed
         type = numeric atom type or type label
         Nsubset = exact number of selected atoms to set to new atom type
         seed = random # seed (positive integer)
       *mol* value = molecule ID
       value can be an atom-style variable (see below)
       *x*,\ *y*,\ *z* value = atom coordinate (distance units)
         value can be an atom-style variable (see below)
       *vx*,\ *vy*,\ *vz* value = atom velocity (velocity units)
         value can be an atom-style variable (see below)
       *charge* value = atomic charge (charge units)
         value can be an atom-style variable (see below)
       *dipole* values = x y z
         x,y,z = orientation of dipole moment vector
         any of x,y,z can be an atom-style variable (see below)
       *dipole/random* value = seed Dlen
         seed = random # seed (positive integer) for dipole moment orientations
         Dlen = magnitude of dipole moment (dipole units)
       *spin/atom* values = g x y z
         g = magnitude of magnetic spin vector (in Bohr magneton's unit)
         x,y,z = orientation of magnetic spin vector
         any of x,y,z can be an atom-style variable (see below)
       *spin/atom/random* value = seed Dlen
         seed = random # seed (positive integer) for magnetic spin orientations
         Dlen = magnitude of magnetic spin vector (in Bohr magneton's unit)
       *radius/electron* values = eradius
         eradius = electron radius (or fixed-core radius) (distance units)
       *spin/electron* value = espin
         espin = electron spin (+1/-1), 0 = nuclei, 2 = fixed-core, 3 = pseudo-cores (i.e. ECP)
       *quat* values = a b c theta
         a,b,c = unit vector to rotate particle around via right-hand rule
         theta = rotation angle (degrees)
         any of a,b,c,theta can be an atom-style variable (see below)
       *quat/random* value = seed
         seed = random # seed (positive integer) for quaternion orientations
       *diameter* value = diameter of spherical particle (distance units)
         value can be an atom-style variable (see below)
       *shape* value = Sx Sy Sz
         Sx,Sy,Sz = 3 diameters of ellipsoid (distance units)
       *length* value = len
         len = length of line segment (distance units)
         len can be an atom-style variable (see below)
       *tri* value = side
         side = side length of equilateral triangle (distance units)
         side can be an atom-style variable (see below)
       *theta* value = angle (degrees)
         angle = orientation of line segment with respect to x-axis
         angle can be an atom-style variable (see below)
       *theta/random* value = seed
         seed = random # seed (positive integer) for line segment orienations
       *angmom* values = Lx Ly Lz
         Lx,Ly,Lz = components of angular momentum vector (distance-mass-velocity units)
         any of Lx,Ly,Lz can be an atom-style variable (see below)
       *omega* values = Wx Wy Wz
         Wx,Wy,Wz = components of angular velocity vector (radians/time units)
         any of wx,wy,wz can be an atom-style variable (see below)
       *mass* value = per-atom mass (mass units)
         value can be an atom-style variable (see below)
       *density* value = particle density for a sphere or ellipsoid (mass/distance\^3 units), or for a triangle (mass/distance\^2 units) or line (mass/distance units) particle
         value can be an atom-style variable (see below)
       *density/disc* value = particle density for a 2d disc or ellipse (mass/distance\^2 units)
         value can be an atom-style variable (see below)
       *temperature* value = temperature for finite-size particles (temperature units)
         value can be an atom-style variable (see below)
       *volume* value = particle volume for Peridynamic particle (distance\^3 units)
         value can be an atom-style variable (see below)
       *image* nx ny nz
         nx,ny,nz = which periodic image of the simulation box the atom is in
         any of nx,ny,nz can be an atom-style variable (see below)
       *bond* value = numeric bond type or bond type label, for all bonds between selected atoms
       *angle* value = numeric angle type or angle type label, for all angles between selected atoms
       *dihedral* value = numeric dihedral type or dihedral type label, for all dihedrals between selected atoms
       *improper* value = numeric improper type or improper type label, for all impropers between selected atoms
       *sph/e* value = energy of SPH particles (need units)
         value can be an atom-style variable (see below)
       *sph/cv* value = heat capacity of SPH particles (need units)
         value can be an atom-style variable (see below)
       *sph/rho* value = density of SPH particles (need units)
         value can be an atom-style variable (see below)
       *smd/contact/radius* = radius for short range interactions, i.e. contact and friction
         value can be an atom-style variable (see below)
       *smd/mass/density* = set particle mass based on volume by providing a mass density
         value can be an atom-style variable (see below)
       *dpd/theta* value = internal temperature of DPD particles (temperature units)
         value can be an atom-style variable (see below)
         value can be NULL which sets internal temp of each particle to KE temp
       *edpd/temp* value = temperature of eDPD particles (temperature units)
         value can be an atom-style variable (see below)
       *edpd/cv* value = volumetric heat capacity of eDPD particles (energy/temperature/volume units)
         value can be an atom-style variable (see below)
       *cc* values = index cc
         index = index of a chemical species (1 to Nspecies)
         cc = chemical concentration of tDPD particles for a species (mole/volume units)
       *epsilon* value = dielectric constant of the medium where the atoms reside
       *i_name* value = custom integer vector with name
       *d_name* value = custom floating-point vector with name
       *i2_name* value = column of a custom integer array with name
                         column specified as i2_name[N] where N is 1 to Ncol
       *d2_name* value = column of a custom floating-point array with name
                         column specified as d2_name[N] where N is 1 to Ncol

Examples
""""""""

.. code-block:: LAMMPS

   set group solvent type 2
   set group solvent type C
   set group solvent type/fraction 2 0.5 12393
   set group solvent type/fraction C 0.5 12393
   set group edge bond 4
   set region half charge 0.5
   set type 3 charge 0.5
   set type H charge 0.5
   set type 1*3 charge 0.5
   set atom * charge v_atomfile
   set atom 100*200 x 0.5 y 1.0
   set atom 100 vx 0.0 vy 0.0 vz -1.0
   set atom 1492 type 3
   set atom 1492 type H
   set atom * i_myVal 5
   set atom * d2_Sxyz[1] 6.4

Description
"""""""""""

Set one or more properties of one or more atoms.  Since atom
properties are initially assigned by the :doc:`read_data <read_data>`,
:doc:`read_restart <read_restart>` or :doc:`create_atoms <create_atoms>`
commands, this command changes those assignments.  This can be useful
for overriding the default values assigned by the
:doc:`create_atoms <create_atoms>` command (e.g. charge = 0.0).  It can
be useful for altering pairwise and molecular force interactions,
since force-field coefficients are defined in terms of types.  It can
be used to change the labeling of atoms by atom type or molecule ID
when they are output in :doc:`dump <dump>` files.  It can also be useful
for debugging purposes; i.e. positioning an atom at a precise location
to compute subsequent forces or energy.

Note that the *style* and *ID* arguments determine which atoms have
their properties reset.  The remaining keywords specify which
properties to reset and what the new values are.  Some strings like
*type* or *mol* can be used as a style and/or a keyword.

----------

This section describes how to select which atoms to change
the properties of, via the *style* and *ID* arguments.

.. versionchanged:: 28Mar2023

   Support for type labels was added for selecting atoms by type

The style *atom* selects all the atoms in a range of atom IDs.

The style *type* selects all the atoms in a range of types or type
labels.  The style *type* selects atoms in one of two ways.  A range
of numeric atom types can be specified.  Or a single atom type label
can be specified, e.g. "C".  The style *mol* selects all the atoms in
a range of molecule IDs.

In each of the range cases, the range can be specified as a single
numeric value, or a wildcard asterisk can be used to specify a range
of values.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  For
example, for the style *type*, if N = the number of atom types, then
an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).  For all the styles except
*mol*, the lowest value for the wildcard is 1; for *mol* it is 0.

The style *group* selects all the atoms in the specified group.  The
style *region* selects all the atoms in the specified geometric
region.  See the :doc:`group <group>` and :doc:`region <region>` commands
for details of how to specify a group or region.

----------

This section describes the keyword options for which properties to
change, for the selected atoms.

Note that except where explicitly prohibited below, all of the
keywords allow an :doc:`atom-style or atomfile-style variable
<variable>` to be used as the specified value(s).  If the value is a
variable, it should be specified as v_name, where name is the
variable name.  In this case, the variable will be evaluated, and its
resulting per-atom value used to determine the value assigned to each
selected atom.  Note that the per-atom value from the variable will be
ignored for atoms that are not selected via the *style* and *ID*
settings explained above.  A simple way to use per-atom values from
the variable to reset a property for all atoms is to use style *atom*
with *ID* = "\*"; this selects all atom IDs.

Atom-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  They can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a time-dependent or
spatially-dependent set of per-atom values.  As explained on the
:doc:`variable <variable>` doc page, atomfile-style variables can be
used in place of atom-style variables, and thus as arguments to the
set command.  Atomfile-style variables read their per-atoms values
from a file.

.. note::

   Atom-style and atomfile-style variables return floating point
   per-atom values.  If the values are assigned to an integer variable,
   such as the molecule ID, then the floating point value is truncated to
   its integer portion, e.g. a value of 2.6 would become 2.

.. versionchanged:: 28Mar2023

   Support for type labels was added for setting atom, bond, angle,
   dihedral, and improper types

Keyword *type* sets the atom type for all selected atoms.  A specified
value can be either a numeric atom type or an atom type label. When
using a numeric type, the specified value must be from 1 to ntypes,
where ntypes was set by the :doc:`create_box <create_box>` command or
the *atom types* field in the header of the data file read by the
:doc:`read_data <read_data>` command.  When using a type label it must
have been defined previously.  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for the allowed syntax of type labels
and a general discussion of how type labels can be used.

Keyword *type/fraction* sets the atom type for a fraction of the selected
atoms.  The actual number of atoms changed is not guaranteed
to be exactly the specified fraction (0 <= *fraction* <= 1), but
should be statistically close.  Random numbers are used in such a way
that a particular atom is changed or not changed, regardless of how
many processors are being used.  This keyword does not allow use of an
atom-style variable.

Keywords *type/ratio* and *type/subset* also set the atom type for a
fraction of the selected atoms.  The actual number of atoms changed
will be exactly the requested number.  For *type/ratio* the specified
fraction (0 <= *fraction* <= 1) determines the number.  For
*type/subset*, the specified *Nsubset* is the number.  An iterative
algorithm is used which ensures the correct number of atoms are
selected, in a perfectly random fashion.  Which atoms are selected
will change with the number of processors used.  These keywords do not
allow use of an atom-style variable.

Keyword *mol* sets the molecule ID for all selected atoms.  The
:doc:`atom style <atom_style>` being used must support the use of
molecule IDs.

Keywords *x*, *y*, *z*, and *charge* set the coordinates or
charge of all selected atoms.  For *charge*, the :doc:`atom style
<atom_style>` being used must support the use of atomic
charge. Keywords *vx*, *vy*, and *vz* set the velocities of all
selected atoms.

Keyword *dipole* uses the specified x,y,z values as components of a
vector to set as the orientation of the dipole moment vectors of the
selected atoms.  The magnitude of the dipole moment is set by the
length of this orientation vector.

Keyword *dipole/random* randomizes the orientation of the dipole
moment vectors for the selected atoms and sets the magnitude of each
to the specified *Dlen* value.  For 2d systems, the z component of the
orientation is set to 0.0.  Random numbers are used in such a way that
the orientation of a particular atom is the same, regardless of how
many processors are being used.  This keyword does not allow use of an
atom-style variable.

.. versionchanged:: 15Sep2022

Keyword *spin/atom* uses the specified g value to set the magnitude of the
magnetic spin vectors, and the x,y,z values as components of a vector
to set as the orientation of the magnetic spin vectors of the selected
atoms.  This keyword was previously called *spin*.

.. versionchanged:: 15Sep2022

Keyword *spin/atom/random* randomizes the orientation of the magnetic spin
vectors for the selected atoms and sets the magnitude of each to the
specified *Dlen* value.  This keyword was previously called *spin/random*.

.. versionadded:: 15Sep2022

Keyword *radius/electron* uses the specified value to set the radius of
electrons or fixed cores.

.. versionadded:: 15Sep2022

Keyword *spin/electron* sets the spin of an electron (+/- 1) or indicates
nuclei (=0), fixed-cores (=2), or pseudo-cores (= 3).

Keyword *quat* uses the specified values to create a quaternion
(4-vector) that represents the orientation of the selected atoms.  The
particles must define a quaternion for their orientation
(e.g. ellipsoids, triangles, body particles) as defined by the
:doc:`atom_style <atom_style>` command.  Note that particles defined by
:doc:`atom_style ellipsoid <atom_style>` have 3 shape parameters.  The 3
values must be non-zero for each particle set by this command.  They
are used to specify the aspect ratios of an ellipsoidal particle,
which is oriented by default with its x-axis along the simulation
box's x-axis, and similarly for y and z.  If this body is rotated (via
the right-hand rule) by an angle theta around a unit rotation vector
(a,b,c), then the quaternion that represents its new orientation is
given by (cos(theta/2), a\*sin(theta/2), b\*sin(theta/2),
c\*sin(theta/2)).  The theta and a,b,c values are the arguments to the
*quat* keyword.  LAMMPS normalizes the quaternion in case (a,b,c) was
not specified as a unit vector.  For 2d systems, the a,b,c values are
ignored, since a rotation vector of (0,0,1) is the only valid choice.

Keyword *quat/random* randomizes the orientation of the quaternion for
the selected atoms.  The particles must define a quaternion for their
orientation (e.g. ellipsoids, triangles, body particles) as defined by
the :doc:`atom_style <atom_style>` command.  Random numbers are used in
such a way that the orientation of a particular atom is the same,
regardless of how many processors are being used.  For 2d systems,
only orientations in the xy plane are generated.  As with keyword
*quat*, for ellipsoidal particles, the 3 shape values must be non-zero
for each particle set by this command.  This keyword does not allow
use of an atom-style variable.

Keyword *diameter* sets the size of the selected atoms.  The particles
must be finite-size spheres as defined by the :doc:`atom_style sphere
<atom_style>` command.  The diameter of a particle can be set to 0.0,
which means they will be treated as point particles.  Note that this
command does not adjust the particle mass, even if it was defined with
a density, e.g. via the :doc:`read_data <read_data>` command.

Keyword *shape* sets the size and shape of the selected atoms.  The
particles must be ellipsoids as defined by the :doc:`atom_style
ellipsoid <atom_style>` command.  The *Sx*, *Sy*, *Sz* settings
are the 3 diameters of the ellipsoid in each direction.  All 3 can be
set to the same value, which means the ellipsoid is effectively a
sphere.  They can also all be set to 0.0 which means the particle will
be treated as a point particle.  Note that this command does not
adjust the particle mass, even if it was defined with a density,
e.g. via the :doc:`read_data <read_data>` command.

Keyword *length* sets the length of selected atoms.  The particles
must be line segments as defined by the :doc:`atom_style line
<atom_style>` command.  If the specified value is non-zero the line
segment is (re)set to a length = the specified value, centered around
the particle position, with an orientation along the x-axis.  If the
specified value is 0.0, the particle will become a point particle.
Note that this command does not adjust the particle mass, even if it
was defined with a density, e.g. via the :doc:`read_data <read_data>`
command.

Keyword *tri* sets the size of selected atoms.  The particles must be
triangles as defined by the :doc:`atom_style tri <atom_style>` command.
If the specified value is non-zero the triangle is (re)set to be an
equilateral triangle in the xy plane with side length = the specified
value, with a centroid at the particle position, with its base
parallel to the x axis, and the y-axis running from the center of the
base to the top point of the triangle.  If the specified value is 0.0,
the particle will become a point particle.  Note that this command
does not adjust the particle mass, even if it was defined with a
density, e.g. via the :doc:`read_data <read_data>` command.

Keyword *theta* sets the orientation of selected atoms.  The particles
must be line segments as defined by the :doc:`atom_style line
<atom_style>` command.  The specified value is used to set the
orientation angle of the line segments with respect to the x axis.

Keyword *theta/random* randomizes the orientation of theta for the
selected atoms.  The particles must be line segments as defined by the
:doc:`atom_style line <atom_style>` command.  Random numbers are used in
such a way that the orientation of a particular atom is the same,
regardless of how many processors are being used.  This keyword does
not allow use of an atom-style variable.

Keyword *angmom* sets the angular momentum of selected atoms.  The
particles must be ellipsoids as defined by the :doc:`atom_style
ellipsoid <atom_style>` command or triangles as defined by the
:doc:`atom_style tri <atom_style>` command.  The angular momentum
vector of the particles is set to the 3 specified components.

Keyword *omega* sets the angular velocity of selected atoms.  The
particles must be spheres as defined by the :doc:`atom_style sphere
<atom_style>` command.  The angular velocity vector of the particles is
set to the 3 specified components.

Keyword *mass* sets the mass of all selected particles.  The particles
must have a per-atom mass attribute, as defined by the :doc:`atom_style
<atom_style>` command.  See the "mass" command for how to set mass
values on a per-type basis.

Keyword *density* or *density/disc* also sets the mass of all selected
particles, but in a different way.  The particles must have a per-atom
mass attribute, as defined by the :doc:`atom_style <atom_style>`
command.  If the atom has a radius attribute (see :doc:`atom_style
sphere <atom_style>`) and its radius is non-zero, its mass is set from
the density and particle volume for 3d systems (the input density is
assumed to be in mass/distance\^3 units).  For 2d, the default is for
LAMMPS to model particles with a radius attribute as spheres.  However,
if the *density/disc* keyword is used, then they can be modeled as 2d
discs (circles).  Their mass is set from the density and particle area
(the input density is assumed to be in mass/distance\^2 units).

If the atom has a shape attribute (see :doc:`atom_style ellipsoid
<atom_style>`) and its 3 shape parameters are non-zero, then its mass is
set from the density and particle volume (the input density is assumed
to be in mass/distance\^3 units).  The *density/disc* keyword has no
effect; it does not (yet) treat 3d ellipsoids as 2d ellipses.

If the atom has a length attribute (see :doc:`atom_style line
<atom_style>`) and its length is non-zero, then its mass is set from the
density and line segment length (the input density is assumed to be in
mass/distance units).  If the atom has an area attribute (see
:doc:`atom_style tri <atom_style>`) and its area is non-zero, then its
mass is set from the density and triangle area (the input density is
assumed to be in mass/distance\^2 units).

If none of these cases are valid, then the mass is set to the density
value directly (the input density is assumed to be in mass units).

Keyword *temperature* sets the temperature of a finite-size particle.
Currently, only the GRANULAR package supports this attribute. The
temperature must be added using an instance of
:doc:`fix property/atom <fix_property_atom>` The values for the
temperature must be positive.

Keyword *volume* sets the volume of all selected particles.  Currently,
only the :doc:`atom_style peri <atom_style>` command defines particles
with a volume attribute.  Note that this command does not adjust the
particle mass.

Keyword *image* sets which image of the simulation box the atom is
considered to be in.  An image of 0 means it is inside the box as
defined.  A value of 2 means add 2 box lengths to get the true value.  A
value of -1 means subtract 1 box length to get the true value.  LAMMPS
updates these flags as atoms cross periodic boundaries during the
simulation.  The flags can be output with atom snapshots via the
:doc:`dump <dump>` command.  If a value of NULL is specified for any of
nx,ny,nz, then the current image value for that dimension is unchanged.
For non-periodic dimensions only a value of 0 can be specified.  This
command can be useful after a system has been equilibrated and atoms
have diffused one or more box lengths in various directions.  This
command can then reset the image values for atoms so that they are
effectively inside the simulation box, e.g if a diffusion coefficient is
about to be measured via the :doc:`compute msd <compute_msd>` command.
Care should be taken not to reset the image flags of two atoms in a bond
to the same value if the bond straddles a periodic boundary (rather they
should be different by +/- 1).  This will not affect the dynamics of a
simulation, but may mess up analysis of the trajectories if a LAMMPS
diagnostic or your own analysis relies on the image flags to unwrap a
molecule which straddles the periodic box.

Keywords *bond*, *angle*, *dihedral*, and *improper*, set the bond
type (angle type, etc) of all bonds (angles, etc) of selected atoms to
the specified value.  The value can be a numeric type from 1 to
nbondtypes (nangletypes, etc).  Or it can be a type label (bond type
label, angle type label, etc).  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for the allowed syntax of type labels
and a general discussion of how type labels can be used.  All atoms in
a particular bond (angle, etc) must be selected atoms in order for the
change to be made.  The value of nbondtypes (nangletypes, etc) was set
by the *bond types* (\ *angle types*, etc) field in the header of the
data file read by the :doc:`read_data <read_data>` command.  These
keywords do not allow use of an atom-style variable.

Keywords *sph/e*, *sph/cv*, and *sph/rho* set the energy, heat capacity,
and density of smoothed particle hydrodynamics (SPH) particles.  See
`this PDF guide <PDF/SPH_LAMMPS_userguide.pdf>`_ to using SPH in LAMMPS.

Keyword *smd/mass/density* sets the mass of all selected particles, but
it is only applicable to the Smooth Mach Dynamics package MACHDYN.  It
assumes that the particle volume has already been correctly set and
calculates particle mass from the provided mass density value.

Keyword *smd/contact/radius* only applies to simulations with the Smooth
Mach Dynamics package MACHDYN.  Itsets an interaction radius for
computing short-range interactions, e.g. repulsive forces to prevent
different individual physical bodies from penetrating each other. Note
that the SPH smoothing kernel diameter used for computing long range,
nonlocal interactions, is set using the *diameter* keyword.

Keyword *dpd/theta* sets the internal temperature of a DPD particle as
defined by the DPD-REACT package.  If the specified value is a number it
must be >= 0.0.  If the specified value is NULL, then the kinetic
temperature Tkin of each particle is computed as 3/2 k Tkin = KE = 1/2 m
v\^2 = 1/2 m (vx\*vx+vy\*vy+vz\*vz).  Each particle's internal
temperature is set to Tkin.  If the specified value is an atom-style
variable, then the variable is evaluated for each particle.  If a value
>= 0.0, the internal temperature is set to that value.  If it is < 0.0,
the computation of Tkin is performed and the internal temperature is set
to that value.

Keywords *edpd/temp* and *edpd/cv* set the temperature and volumetric
heat capacity of an eDPD particle as defined by the DPD-MESO package.
Currently, only :doc:`atom_style edpd <atom_style>` defines particles
with these attributes. The values for the temperature and heat capacity
must be positive.

Keyword *cc* sets the chemical concentration of a tDPD particle for a
specified species as defined by the DPD-MESO package. Currently, only
:doc:`atom_style tdpd <atom_style>` defines particles with this
attribute. An integer for "index" selects a chemical species (1 to
Nspecies) where Nspecies is set by the atom_style command. The value for
the chemical concentration must be >= 0.0.

Keyword *epsilon* sets the dielectric constant of a particle, precisely
of the medium where the particle resides as defined by the DIELECTRIC
package. Currently, only :doc:`atom_style dielectric <atom_style>`
defines particles with this attribute. The value for the dielectric
constant must be >= 0.0.  Note that the set command with this keyword
will rescale the particle charge accordingly so that the real charge
(e.g., as read from a data file) stays intact. To change the real
charges, one needs to use the set command with the *charge*
keyword. Care must be taken to ensure that the real and scaled charges,
and dielectric constants are consistent.

Keywords *i_name*, *d_name*, *i2_name*, *d2_name* refer to custom
per-atom integer and floating-point vectors or arrays that have been
added via the :doc:`fix property/atom <fix_property_atom>` command.
When that command is used specific names are given to each attribute
which are the "name" portion of these keywords.  For arrays *i2_name*
and *d2_name*, the column of the array must also be included following
the name in brackets: e.g. d2_xyz[2], i2_mySpin[3].

Restrictions
""""""""""""

You cannot set an atom attribute (e.g. *mol* or *q* or *volume*\ ) if
the :doc:`atom_style <atom_style>` does not have that attribute.

This command requires inter-processor communication to coordinate the
setting of bond types (angle types, etc).  This means that your system
must be ready to perform a simulation before using one of these
keywords (force fields set, atom mass set, etc).  This is not
necessary for other keywords.

Using the *region* style with the bond (angle, etc) keywords can give
unpredictable results if there are bonds (angles, etc) that straddle
periodic boundaries.  This is because the region may only extend up to
the boundary and partner atoms in the bond (angle, etc) may have
coordinates outside the simulation box if they are ghost atoms.

Related commands
""""""""""""""""

:doc:`create_box <create_box>`, :doc:`create_atoms <create_atoms>`,
:doc:`read_data <read_data>`

Default
"""""""

none
