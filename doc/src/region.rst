.. index:: region

region command
==============

Syntax
""""""


.. parsed-literal::

   region ID style args keyword arg ...

* ID = user-assigned name for the region
* style = *delete* or *block* or *cone* or *cylinder* or *plane* or *prism* or *sphere* or *union* or *intersect*

  .. parsed-literal::

       *delete* = no args
       *block* args = xlo xhi ylo yhi zlo zhi
         xlo,xhi,ylo,yhi,zlo,zhi = bounds of block in all dimensions (distance units)
       *cone* args = dim c1 c2 radlo radhi lo hi
         dim = *x* or *y* or *z* = axis of cone
         c1,c2 = coords of cone axis in other 2 dimensions (distance units)
         radlo,radhi = cone radii at lo and hi end (distance units)
         lo,hi = bounds of cone in dim (distance units)
       *cylinder* args = dim c1 c2 radius lo hi
         dim = *x* or *y* or *z* = axis of cylinder
         c1,c2 = coords of cylinder axis in other 2 dimensions (distance units)
         radius = cylinder radius (distance units)
           c1,c2, and radius can be a variable (see below)
         lo,hi = bounds of cylinder in dim (distance units)
       *plane* args = px py pz nx ny nz
         px,py,pz = point on the plane (distance units)
         nx,ny,nz = direction normal to plane (distance units)
       *prism* args = xlo xhi ylo yhi zlo zhi xy xz yz
         xlo,xhi,ylo,yhi,zlo,zhi = bounds of untilted prism (distance units)
         xy = distance to tilt y in x direction (distance units)
         xz = distance to tilt z in x direction (distance units)
         yz = distance to tilt z in y direction (distance units)
       *sphere* args = x y z radius
         x,y,z = center of sphere (distance units)
         radius = radius of sphere (distance units)
           x,y,z, and radius can be a variable (see below)
       *union* args = N reg-ID1 reg-ID2 ...
         N = # of regions to follow, must be 2 or greater
         reg-ID1,reg-ID2, ... = IDs of regions to join together
       *intersect* args = N reg-ID1 reg-ID2 ...
         N = # of regions to follow, must be 2 or greater
         reg-ID1,reg-ID2, ... = IDs of regions to intersect

* zero or more keyword/arg pairs may be appended
* keyword = *side* or *units* or *move* or *rotate* or *open*

  .. parsed-literal::

       *side* value = *in* or *out*
         *in* = the region is inside the specified geometry
         *out* = the region is outside the specified geometry
       *units* value = *lattice* or *box*
         *lattice* = the geometry is defined in lattice units
         *box* = the geometry is defined in simulation box units
       *move* args = v_x v_y v_z
         v_x,v_y,v_z = equal-style variables for x,y,z displacement of region over time
       *rotate* args = v_theta Px Py Pz Rx Ry Rz
         v_theta = equal-style variable for rotaton of region over time (in radians)
         Px,Py,Pz = origin for axis of rotation (distance units)
         Rx,Ry,Rz = axis of rotation vector
       *open* value = integer from 1-6 corresponding to face index (see below)

* accelerated styles (with same args) = *block/kk*


Examples
""""""""


.. parsed-literal::

   region 1 block -3.0 5.0 INF 10.0 INF INF
   region 2 sphere 0.0 0.0 0.0 5 side out
   region void cylinder y 2 3 5 -5.0 EDGE units box
   region 1 prism 0 10 0 10 0 10 2 0 0
   region outside union 4 side1 side2 side3 side4
   region 2 sphere 0.0 0.0 0.0 5 side out move v_left v_up NULL
   region openbox block 0 10 0 10 0 10 open 5 open 6 units box
   region funnel cone z 10 10 2 5 0 10 open 1 units box

Description
"""""""""""

This command defines a geometric region of space.  Various other
commands use regions.  For example, the region can be filled with
atoms via the :doc:`create_atoms <create_atoms>` command.  Or a bounding
box around the region, can be used to define the simulation box via
the :doc:`create_box <create_box>` command.  Or the atoms in the region
can be identified as a group via the :doc:`group <group>` command, or
deleted via the :doc:`delete_atoms <delete_atoms>` command.  Or the
surface of the region can be used as a boundary wall via the :doc:`fix wall/region <fix_wall_region>` command.

Commands which use regions typically test whether an atom's position
is contained in the region or not.  For this purpose, coordinates
exactly on the region boundary are considered to be interior to the
region.  This means, for example, for a spherical region, an atom on
the sphere surface would be part of the region if the sphere were
defined with the *side in* keyword, but would not be part of the
region if it were defined using the *side out* keyword.  See more
details on the *side* keyword below.

Normally, regions in LAMMPS are "static", meaning their geometric
extent does not change with time.  If the *move* or *rotate* keyword
is used, as described below, the region becomes "dynamic", meaning
it's location or orientation changes with time.  This may be useful,
for example, when thermostatting a region, via the compute temp/region
command, or when the fix wall/region command uses a region surface as
a bounding wall on particle motion, i.e. a rotating container.

The *delete* style removes the named region.  Since there is little
overhead to defining extra regions, there is normally no need to do
this, unless you are defining and discarding large numbers of regions
in your input script.

The lo/hi values for *block* or *cone* or *cylinder* or *prism* styles
can be specified as EDGE or INF.  EDGE means they extend all the way
to the global simulation box boundary.  Note that this is the current
box boundary; if the box changes size during a simulation, the region
does not.  INF means a large negative or positive number (1.0e20), so
it should encompass the simulation box even if it changes size.  If a
region is defined before the simulation box has been created (via
:doc:`create_box <create_box>` or :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands), then an EDGE or INF
parameter cannot be used.  For a *prism* region, a non-zero tilt
factor in any pair of dimensions cannot be used if both the lo/hi
values in either of those dimensions are INF.  E.g. if the xy tilt is
non-zero, then xlo and xhi cannot both be INF, nor can ylo and yhi.

.. note::

   Regions in LAMMPS do not get wrapped across periodic boundaries,
   as specified by the :doc:`boundary <boundary>` command.  For example, a
   spherical region that is defined so that it overlaps a periodic
   boundary is not treated as 2 half-spheres, one on either side of the
   simulation box.

.. note::

   Regions in LAMMPS are always 3d geometric objects, regardless of
   whether the :doc:`dimension <dimension>` of a simulation is 2d or 3d.
   Thus when using regions in a 2d simulation, you should be careful to
   define the region so that its intersection with the 2d x-y plane of
   the simulation has the 2d geometric extent you want.

For style *cone*\ , an axis-aligned cone is defined which is like a
*cylinder* except that two different radii (one at each end) can be
defined.  Either of the radii (but not both) can be 0.0.

For style *cone* and *cylinder*\ , the c1,c2 params are coordinates in
the 2 other dimensions besides the cylinder axis dimension.  For dim =
x, c1/c2 = y/z; for dim = y, c1/c2 = x/z; for dim = z, c1/c2 = x/y.
Thus the third example above specifies a cylinder with its axis in the
y-direction located at x = 2.0 and z = 3.0, with a radius of 5.0, and
extending in the y-direction from -5.0 to the upper box boundary.

For style *plane*\ , a plane is defined which contain the point
(px,py,pz) and has a normal vector (nx,ny,nz).  The normal vector does
not have to be of unit length.  The "inside" of the plane is the
half-space in the direction of the normal vector; see the discussion
of the *side* option below.

For style *prism*\ , a parallelepiped is defined (it's too hard to spell
parallelepiped in an input script!).  The parallelepiped has its
"origin" at (xlo,ylo,zlo) and is defined by 3 edge vectors starting
from the origin given by A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C =
(xz,yz,zhi-zlo).  *Xy,xz,yz* can be 0.0 or positive or negative values
and are called "tilt factors" because they are the amount of
displacement applied to faces of an originally orthogonal box to
transform it into the parallelepiped.

A prism region that will be used with the :doc:`create_box <create_box>`
command to define a triclinic simulation box must have tilt factors
(xy,xz,yz) that do not skew the box more than half the distance of
corresponding the parallel box length.  For example, if xlo = 2 and
xhi = 12, then the x box length is 10 and the xy tilt factor must be
between -5 and 5.  Similarly, both xz and yz must be between
-(xhi-xlo)/2 and +(yhi-ylo)/2.  Note that this is not a limitation,
since if the maximum tilt factor is 5 (as in this example), then
configurations with tilt = ..., -15, -5, 5, 15, 25, ... are all
geometrically equivalent.

The *radius* value for style *sphere* and *cylinder* can be specified
as an equal-style :doc:`variable <variable>`.  If the value is a
variable, it should be specified as v\_name, where name is the variable
name.  In this case, the variable will be evaluated each timestep, and
its value used to determine the radius of the region. For style *sphere*
also the x-, y-, and z- coordinate of the center of the sphere and for
style *cylinder* the two center positions c1 and c2 for the location of
the cylinder axes can be a variable with the same kind of effect and
requirements than for the radius.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent radius or have
a time dependent position of the sphere or cylinder region.

See the :doc:`Howto tricilinc <Howto_triclinic>` doc page for a
geometric description of triclinic boxes, as defined by LAMMPS, and
how to transform these parameters to and from other commonly used
triclinic representations.

The *union* style creates a region consisting of the volume of all the
listed regions combined.  The *intersect* style creates a region
consisting of the volume that is common to all the listed regions.

.. note::

   The *union* and *intersect* regions operate by invoking methods
   from their list of sub-regions.  Thus you cannot delete the
   sub-regions after defining a *union* or *intersection* region.


----------


The *side* keyword determines whether the region is considered to be
inside or outside of the specified geometry.  Using this keyword in
conjunction with *union* and *intersect* regions, complex geometries
can be built up.  For example, if the interior of two spheres were
each defined as regions, and a *union* style with *side* = out was
constructed listing the region-IDs of the 2 spheres, the resulting
region would be all the volume in the simulation box that was outside
both of the spheres.

The *units* keyword determines the meaning of the distance units used
to define the region for any argument above listed as having distance
units.  It also affects the scaling of the velocity vector specified
with the *vel* keyword, the amplitude vector specified with the
*wiggle* keyword, and the rotation point specified with the *rotate*
keyword, since they each involve a distance metric.

A *box* value selects standard distance units as defined by the
:doc:`units <units>` command, e.g. Angstroms for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacings which are used as follows:

* For style *block*\ , the lattice spacing in dimension x is applied to
  xlo and xhi, similarly the spacings in dimensions y,z are applied to
  ylo/yhi and zlo/zhi.
* For style *cone*\ , the lattice spacing in argument *dim* is applied to
  lo and hi.  The spacings in the two radial dimensions are applied to
  c1 and c2.  The two cone radii are scaled by the lattice
  spacing in the dimension corresponding to c1.
* For style *cylinder*\ , the lattice spacing in argument *dim* is applied
  to lo and hi.  The spacings in the two radial dimensions are applied
  to c1 and c2.  The cylinder radius is scaled by the lattice
  spacing in the dimension corresponding to c1.
* For style *plane*\ , the lattice spacing in dimension x is applied to
  px and nx, similarly the spacings in dimensions y,z are applied to
  py/ny and pz/nz.
* For style *prism*\ , the lattice spacing in dimension x is applied to
  xlo and xhi, similarly for ylo/yhi and zlo/zhi.  The lattice spacing
  in dimension x is applied to xy and xz, and the spacing in dimension y
  to yz.
* For style *sphere*\ , the lattice spacing in dimensions x,y,z are
  applied to the sphere center x,y,z.  The spacing in dimension x is
  applied to the sphere radius.


----------


If the *move* or *rotate* keywords are used, the region is "dynamic",
meaning its location or orientation changes with time.  These keywords
cannot be used with a *union* or *intersect* style region.  Instead,
the keywords should be used to make the individual sub-regions of the
*union* or *intersect* region dynamic.  Normally, each sub-region
should be "dynamic" in the same manner (e.g. rotate around the same
point), though this is not a requirement.

The *move* keyword allows one or more :doc:`equal-style variables <variable>` to be used to specify the x,y,z displacement
of the region, typically as a function of time.  A variable is
specified as v\_name, where name is the variable name.  Any of the
three variables can be specified as NULL, in which case no
displacement is calculated in that dimension.

Note that equal-style variables can specify formulas with various
mathematical functions, and include :doc:`thermo_style <thermo_style>`
command keywords for the simulation box parameters and timestep and
elapsed time.  Thus it is easy to specify a region displacement that
change as a function of time or spans consecutive runs in a continuous
fashion.  For the latter, see the *start* and *stop* keywords of the
:doc:`run <run>` command and the *elaplong* keyword of :doc:`thermo_style custom <thermo_style>` for details.

For example, these commands would displace a region from its initial
position, in the positive x direction, effectively at a constant
velocity:


.. parsed-literal::

   variable dx equal ramp(0,10)
   region 2 sphere 10.0 10.0 0.0 5 move v_dx NULL NULL

Note that the initial displacement is 0.0, though that is not required.

Either of these variables would "wiggle" the region back and forth in
the y direction:


.. parsed-literal::

   variable dy equal swiggle(0,5,100)
   variable dysame equal 5\*sin(2\*PI\*elaplong\*dt/100)
   region 2 sphere 10.0 10.0 0.0 5 move NULL v_dy NULL

The *rotate* keyword rotates the region around a rotation axis *R* =
(Rx,Ry,Rz) that goes through a point *P* = (Px,Py,Pz).  The rotation
angle is calculated, presumably as a function of time, by a variable
specified as v\_theta, where theta is the variable name.  The variable
should generate its result in radians.  The direction of rotation for
the region around the rotation axis is consistent with the right-hand
rule: if your right-hand thumb points along *R*\ , then your fingers
wrap around the axis in the direction of rotation.

The *move* and *rotate* keywords can be used together.  In this case,
the displacement specified by the *move* keyword is applied to the *P*
point of the *rotate* keyword.


----------


The *open* keyword can be used (multiple times) to indicate that one
or more faces of the region are ignored for purposes of particle/wall
interactions.  This keyword is only relevant for regions used by the
*fix wall/region* and *fix wall/gran/region* commands.  It can be used
to create "open" containers where only some of the region faces are
walls.  For example, a funnel can be created with a *cone* style
region that has an open face at the smaller radius for particles to
flow out, or at the larger radius for pouring particles into the cone,
or both.

Note that using the *open* keyword partly overrides the *side*
keyword, since both exterior and interior surfaces of an open region
are tested for particle contacts.  The exception to this is a *union*
or *intersect* region which includes an open sub-region.  In that case
the *side* keyword is still used to define the union/intersect region
volume, and the *open* settings are only applied to the individual
sub-regions that use them.

The indices specified as part of the *open* keyword have the following
meanings:

For style *block*\ , indices 1-6 correspond to the xlo, xhi, ylo, yhi,
zlo, zhi surfaces of the block.  I.e. 1 is the yz plane at x = xlo, 2
is the yz-plane at x = xhi, 3 is the xz plane at y = ylo, 4 is the xz
plane at y = yhi, 5 is the xy plane at z = zlo, 6 is the xy plane at z
= zhi).  In the second-to-last example above, the region is a box open
at both xy planes.

For style *prism*\ , values 1-6 have the same mapping as for style
*block*\ .  I.e. in an untilted *prism*\ , *open* indices correspond to
the xlo, xhi, ylo, yhi, zlo, zhi surfaces.

For style *cylinder*\ , index 1 corresponds to the flat end cap at the
low coordinate along the cylinder axis, index 2 corresponds to the
high-coordinate flat end cap along the cylinder axis, and index 3 is
the curved cylinder surface.  For example, a *cylinder* region with
*open 1 open 2* keywords will be open at both ends (e.g. a section of
pipe), regardless of the cylinder orientation.

For style *cone*\ , the mapping is the same as for style *cylinder*\ .
Index 1 is the low-coordinate flat end cap, index 2 is the
high-coordinate flat end cap, and index 3 is the curved cone surface.
In the last example above, a *cone* region is defined along the z-axis
that is open at the zlo value (e.g. for use as a funnel).

For all other styles, the *open* keyword is ignored.  As indicated
above, this includes the *intersect* and *union* regions, though their
sub-regions can be defined with the *open* keyword.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

The code using the region (such as a fix or compute) must also be supported
by Kokkos or no acceleration will occur. Currently, only *block* style
regions are supported by Kokkos.

These accelerated styles are part of the Kokkos package.  They are
only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


A prism cannot be of 0.0 thickness in any dimension; use a small z
thickness for 2d simulations.  For 2d simulations, the xz and yz
parameters must be 0.0.

Related commands
""""""""""""""""

:doc:`lattice <lattice>`, :doc:`create_atoms <create_atoms>`,
:doc:`delete_atoms <delete_atoms>`, :doc:`group <group>`

Default
"""""""

The option defaults are side = in, units = lattice, and no move or
rotation.
