.. index:: lattice

lattice command
===============

Syntax
""""""

.. parsed-literal::

   lattice style scale keyword values ...

* style = *none* or *sc* or *bcc* or *fcc* or *hcp* or *diamond* or         *sq* or *sq2* or *hex* or *custom*
* scale = scale factor between lattice and simulation box

  .. parsed-literal::

       scale = reduced density rho\* (for LJ units)
       scale = lattice constant in distance units (for all other units)

* zero or more keyword/value pairs may be appended
* keyword = *origin* or *orient* or *spacing* or *a1* or *a2* or *a3* or *basis*

  .. parsed-literal::

       *origin* values = x y z
         x,y,z = fractions of a unit cell (0 <= x,y,z < 1)
       *orient* values = dim i j k
         dim = *x* or *y* or *z*
         i,j,k = integer lattice directions
       *spacing* values = dx dy dz
         dx,dy,dz = lattice spacings in the x,y,z box directions
       *a1*\ ,\ *a2*\ ,\ *a3* values = x y z
         x,y,z = primitive vector components that define unit cell
       *basis* values = x y z
         x,y,z = fractional coords of a basis atom (0 <= x,y,z < 1)

Examples
""""""""

.. parsed-literal::

   lattice fcc 3.52
   lattice hex 0.85
   lattice sq 0.8 origin 0.0 0.5 0.0 orient x 1 1 0 orient y -1 1 0
   lattice custom 3.52 a1 1.0 0.0 0.0 a2 0.5 1.0 0.0 a3 0.0 0.0 0.5 &
                       basis 0.0 0.0 0.0 basis 0.5 0.5 0.5
   lattice none 2.0

Description
"""""""""""

Define a lattice for use by other commands.  In LAMMPS, a lattice is
simply a set of points in space, determined by a unit cell with basis
atoms, that is replicated infinitely in all dimensions.  The arguments
of the lattice command can be used to define a wide variety of
crystallographic lattices.

A lattice is used by LAMMPS in two ways.  First, the
:doc:`create_atoms <create_atoms>` command creates atoms on the lattice
points inside the simulation box.  Note that the
:doc:`create_atoms <create_atoms>` command allows different atom types
to be assigned to different basis atoms of the lattice.  Second, the
lattice spacing in the x,y,z dimensions implied by the lattice, can be
used by other commands as distance units
(e.g. :doc:`create_box <create_box>`, :doc:`region <region>` and
:doc:`velocity <velocity>`), which are often convenient to use when the
underlying problem geometry is atoms on a lattice.

The lattice style must be consistent with the dimension of the
simulation - see the :doc:`dimension <dimension>` command.  Styles *sc*
or *bcc* or *fcc* or *hcp* or *diamond* are for 3d problems.  Styles
*sq* or *sq2* or *hex* are for 2d problems.  Style *custom* can be
used for either 2d or 3d problems.

A lattice consists of a unit cell, a set of basis atoms within that
cell, and a set of transformation parameters (scale, origin, orient)
that map the unit cell into the simulation box.  The vectors a1,a2,a3
are the edge vectors of the unit cell.  This is the nomenclature for
"primitive" vectors in solid-state crystallography, but in LAMMPS the
unit cell they determine does not have to be a "primitive cell" of
minimum volume.

Note that the lattice command can be used multiple times in an input
script.  Each time it is invoked, the lattice attributes are
re-defined and are used for all subsequent commands (that use lattice
attributes).  For example, a sequence of lattice,
:doc:`region <region>`, and :doc:`create_atoms <create_atoms>` commands
can be repeated multiple times to build a poly-crystalline model with
different geometric regions populated with atoms in different lattice
orientations.

----------

A lattice of style *none* does not define a unit cell and basis set,
so it cannot be used with the :doc:`create_atoms <create_atoms>`
command.  However it does define a lattice spacing via the specified
scale parameter.  As explained above the lattice spacings in x,y,z can
be used by other commands as distance units.  No additional
keyword/value pairs can be specified for the *none* style.  By
default, a "lattice none 1.0" is defined, which means the lattice
spacing is the same as one distance unit, as defined by the
:doc:`units <units>` command.

Lattices of style *sc*\ , *fcc*\ , *bcc*\ , and *diamond* are 3d lattices
that define a cubic unit cell with edge length = 1.0.  This means a1 =
1 0 0, a2 = 0 1 0, and a3 = 0 0 1.  Style *hcp* has a1 = 1 0 0, a2 = 0
sqrt(3) 0, and a3 = 0 0 sqrt(8/3).  The placement of the basis atoms
within the unit cell are described in any solid-state physics text.  A
*sc* lattice has 1 basis atom at the lower-left-bottom corner of the
cube.  A *bcc* lattice has 2 basis atoms, one at the corner and one at
the center of the cube.  A *fcc* lattice has 4 basis atoms, one at the
corner and 3 at the cube face centers.  A *hcp* lattice has 4 basis
atoms, two in the z = 0 plane and 2 in the z = 0.5 plane.  A *diamond*
lattice has 8 basis atoms.

Lattices of style *sq* and *sq2* are 2d lattices that define a square
unit cell with edge length = 1.0.  This means a1 = 1 0 0 and a2 = 0 1
0.  A *sq* lattice has 1 basis atom at the lower-left corner of the
square.  A *sq2* lattice has 2 basis atoms, one at the corner and one
at the center of the square.  A *hex* style is also a 2d lattice, but
the unit cell is rectangular, with a1 = 1 0 0 and a2 = 0 sqrt(3) 0.
It has 2 basis atoms, one at the corner and one at the center of the
rectangle.

A lattice of style *custom* allows you to specify a1, a2, a3, and a
list of basis atoms to put in the unit cell.  By default, a1 and a2
and a3 are 3 orthogonal unit vectors (edges of a unit cube).  But you
can specify them to be of any length and non-orthogonal to each other,
so that they describe a tilted parallelepiped.  Via the *basis*
keyword you add atoms, one at a time, to the unit cell.  Its arguments
are fractional coordinates (0.0 <= x,y,z < 1.0).  The position vector
x of a basis atom within the unit cell is thus a linear combination of
the unit cell's 3 edge vectors, i.e. x = bx a1 + by a2 + bz a3,
where bx,by,bz are the 3 values specified for the *basis* keyword.

----------

This sub-section discusses the arguments that determine how the
idealized unit cell is transformed into a lattice of points within the
simulation box.

The *scale* argument determines how the size of the unit cell will be
scaled when mapping it into the simulation box.  I.e. it determines a
multiplicative factor to apply to the unit cell, to convert it to a
lattice of the desired size and distance units in the simulation box.
The meaning of the *scale* argument depends on the :doc:`units <units>`
being used in your simulation.

For all unit styles except *lj*\ , the scale argument is specified in
the distance units defined by the unit style.  For example, in *real*
or *metal* units, if the unit cell is a unit cube with edge length
1.0, specifying scale = 3.52 would create a cubic lattice with a
spacing of 3.52 Angstroms.  In *cgs* units, the spacing would be 3.52
cm.

For unit style *lj*\ , the scale argument is the Lennard-Jones reduced
density, typically written as rho\*.  LAMMPS converts this value into
the multiplicative factor via the formula "factor\^dim = rho/rho\*",
where rho = N/V with V = the volume of the lattice unit cell and N =
the number of basis atoms in the unit cell (described below), and dim
= 2 or 3 for the dimensionality of the simulation.  Effectively, this
means that if LJ particles of size sigma = 1.0 are used in the
simulation, the lattice of particles will be at the desired reduced
density.

The *origin* option specifies how the unit cell will be shifted or
translated when mapping it into the simulation box.  The x,y,z values
are fractional values (0.0 <= x,y,z < 1.0) meaning shift the lattice
by a fraction of the lattice spacing in each dimension.  The meaning
of "lattice spacing" is discussed below.

The *orient* option specifies how the unit cell will be rotated when
mapping it into the simulation box.  The *dim* argument is one of the
3 coordinate axes in the simulation box.  The other 3 arguments are
the crystallographic direction in the lattice that you want to orient
along that axis, specified as integers.  E.g. "orient x 2 1 0" means
the x-axis in the simulation box will be the [210] lattice
direction, and similarly for y and z.  The 3 lattice directions you
specify do not have to be unit vectors, but they must be mutually
orthogonal and obey the right-hand rule, i.e. (X cross Y) points in
the Z direction.

.. note::

   The preceding paragraph describing lattice directions is only
   valid for orthogonal cubic unit cells (or square in 2d).  If you are
   using a *hcp* or *hex* lattice or the more general lattice style
   *custom* with non-orthogonal a1,a2,a3 vectors, then you should think
   of the 3 *orient* vectors as creating a 3x3 rotation matrix which is
   applied to a1,a2,a3 to rotate the original unit cell to a new
   orientation in the simulation box.

----------

Several LAMMPS commands have the option to use distance units that are
inferred from "lattice spacings" in the x,y,z box directions.
E.g. the :doc:`region <region>` command can create a block of size
10x20x20, where 10 means 10 lattice spacings in the x direction.

.. note::

   Though they are called lattice spacings, all the commands that
   have a "units lattice" option, simply use the 3 values as scale
   factors on the distance units defined by the :doc:`units <units>`
   command.  Thus if you do not like the lattice spacings computed by
   LAMMPS (e.g. for a non-orthogonal or rotated unit cell), you can
   define the 3 values to be whatever you wish, via the *spacing* option.

If the *spacing* option is not specified, the lattice spacings are
computed by LAMMPS in the following way.  A unit cell of the lattice
is mapped into the simulation box (scaled and rotated), so that it now
has (perhaps) a modified size and orientation.  The lattice spacing in
X is defined as the difference between the min/max extent of the x
coordinates of the 8 corner points of the modified unit cell (4 in
2d).  Similarly, the Y and Z lattice spacings are defined as the
difference in the min/max of the y and z coordinates.

Note that if the unit cell is orthogonal with axis-aligned edges (no
rotation via the *orient* keyword), then the lattice spacings in each
dimension are simply the scale factor (described above) multiplied by
the length of a1,a2,a3.  Thus a *hex* style lattice with a scale
factor of 3.0 Angstroms, would have a lattice spacing of 3.0 in x and
3\*sqrt(3.0) in y.

.. note::

   For non-orthogonal unit cells and/or when a rotation is applied
   via the *orient* keyword, then the lattice spacings computed by LAMMPS
   are typically less intuitive.  In particular, in these cases, there is
   no guarantee that a particular lattice spacing is an integer multiple
   of the periodicity of the lattice in that direction.  Thus, if you
   create an orthogonal periodic simulation box whose size in a dimension
   is a multiple of the lattice spacing, and then fill it with atoms via
   the :doc:`create_atoms <create_atoms>` command, you will NOT necessarily
   create a periodic system.  I.e. atoms may overlap incorrectly at the
   faces of the simulation box.

The *spacing* option sets the 3 lattice spacings directly.  All must
be non-zero (use 1.0 for dz in a 2d simulation).  The specified values
are multiplied by the multiplicative factor described above that is
associated with the scale factor.  Thus a spacing of 1.0 means one
unit cell edge length independent of the scale factor.  As mentioned
above, this option can be useful if the spacings LAMMPS computes are
inconvenient to use in subsequent commands, which can be the case for
non-orthogonal or rotated lattices.

Note that whenever the lattice command is used, the values of the
lattice spacings LAMMPS calculates are printed out.  Thus their effect
in commands that use the spacings should be decipherable.

----------

Example commands for generating a Wurtzite crystal (courtesy
of Aidan Thompson), with its 8 atom unit cell.

.. parsed-literal::

   variable a equal  4.340330
   variable b equal  $a\*sqrt(3.0)
   variable c equal  $a\*sqrt(8.0/3.0)

   variable 1_3 equal 1.0/3.0
   variable 2_3 equal 2.0/3.0
   variable 1_6 equal 1.0/6.0
   variable 5_6 equal 5.0/6.0
   variable 1_12 equal 1.0/12.0
   variable 5_12 equal 5.0/12.0

   lattice custom    1.0     &
           a1      $a      0.0     0.0     &
           a2      0.0     $b      0.0     &
           a3      0.0     0.0     $c      &
           basis   0.0     0.0     0.0     &
           basis   0.5     0.5     0.0     &
           basis   ${1_3}  0.0     0.5     &
           basis   ${5_6}  0.5     0.5     &
           basis   0.0     0.0     0.625   &
           basis   0.5     0.5     0.625   &
           basis   ${1_3}  0.0     0.125   &
           basis   ${5_6}  0.5     0.125

   region myreg block 0 1 0 1 0 1
   create_box      2 myreg
   create_atoms    1 box      &
           basis   5       2       &
           basis   6       2       &
           basis   7       2       &
           basis   8       2

----------

Restrictions
""""""""""""

The *a1,a2,a3,basis* keywords can only be used with style *custom*\ .

Related commands
""""""""""""""""

:doc:`dimension <dimension>`, :doc:`create_atoms <create_atoms>`,
:doc:`region <region>`

Default
"""""""

.. parsed-literal::

   lattice none 1.0

For other lattice styles, the option defaults are origin = 0.0 0.0
0.0, orient = x 1 0 0, orient = y 0 1 0, orient = z 0 0 1, a1 = 1 0 0,
a2 = 0 1 0, and a3 = 0 0 1.
