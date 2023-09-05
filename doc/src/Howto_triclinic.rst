Triclinic (non-orthogonal) simulation boxes
===========================================

By default, LAMMPS uses an orthogonal simulation box to encompass the
particles.  The orthogonal box has its "origin" at (xlo,ylo,zlo) and
extends to (xhi,yhi,zhi).  Conceptually it is defined by 3 edge
vectors starting from the origin given by **A** = (xhi-xlo,0,0); **B**
= (0,yhi-ylo,0); **C** = (0,0,zhi-zlo).  The :doc:`boundary
<boundary>` command sets the boundary conditions for the 6 faces of
the box (periodic, non-periodic, etc).  The 6 parameters
(xlo,xhi,ylo,yhi,zlo,zhi) are defined at the time the simulation box
is created by one of these commands:

* :doc:`create_box <create_box>`
* :doc:`read_data <read_data>`
* :doc:`read_restart <read_restart>`
* :doc:`read_dump <read_dump>`

Internally, LAMMPS defines box size parameters lx,ly,lz where lx =
xhi-xlo, and similarly in the y and z dimensions.  The 6 parameters,
as well as lx,ly,lz, can be output via the :doc:`thermo_style custom
<thermo_style>` command.  See the :doc:'Howto 2d <Howto_2d>` doc page
for info on how zlo and zhi are defined for 2d simulations.

----------

LAMMPS also allows simulations to be performed using triclinic
(non-orthogonal) simulation boxes shaped as a parallelepiped with
triclinic symmetry.

One use of triclinic simulation boxes is to model solid-state crystals
with triclinic symmetry.  The :doc:`lattice <lattice>` command can be
used with non-orthogonal basis vectors to define a lattice that will
tile a triclinic simulation box via the :doc:`create_atoms
<create_atoms>` command.

A second use is to run Parrinello-Rahman dynamics via the :doc:`fix
npt <fix_nh>` command, which will adjust the xy, xz, yz tilt factors
to compensate for off-diagonal components of the pressure tensor.  The
analog for an :doc:`energy minimization <minimize>` is the :doc:`fix
box/relax <fix_box_relax>` command.

A third use is to shear a bulk solid to study the response of the
material.  The :doc:`fix deform <fix_deform>` command can be used for
this purpose.  It allows dynamic control of the xy, xz, yz tilt
factors as a simulation runs.  This is discussed in the next section
on non-equilibrium MD (NEMD) simulations.

Conceptually, the tricliic parallelepiped is defined with an "origin"
at (xlo,ylo,zhi) and 3 edge vectors **A** = (ax,ay,az), **B** =
(bx,by,bz), **C** = (cx,cy,cz) which can now be arbitrary vectors, so
long as they are distinct and are not co-planar.

The 4 commands listed above for defining orthogonal simulation boxes
have triclinic options which allow for specification of the origin and
edge vectors **A**, **B**, **C**.  For each command, this can be done
in one of two ways, for what LAMMPS calls a *general* triclinic box,
or a *restricited* triclinic box.  A *general* triclinic box is
specified by an origin and 9 parameters (ax,ay,az), (bx,by,bz),
(cx,cy,cz), or 12 parameters in total.  A *restricted* triclinic box
also has an origin, but its edge vectors are of the following form:
**A** = (xhi-xlo,0,0), **B** = (xy,yhi-ylo,0), **C** =
(xz,yz,zhi-zlo).  So 9 parameters in total.

The restricted form of edge vectors means that **A** is along the
x-axis, **B** is in the x-y plane with a y-component in the +y
direction, and **C** has a z-component in the +z direction.
*Xy,xz,yz* can be 0.0 or positive or negative values and are called
"tilt factors" because they are the amount of displacement applied to
faces of an originally orthogonal box to transform it into a
restricted parallelepiped.

.. note::

   Any general triclinic box (i.e. any solid-state crystal basis
   vectors) can be rotated/inverted in 3d around its origin to conform
   to the definition of a restricted triclinic box.  An inversion may
   need to be applied to the rotated **C** vector to ensure its final
   z-component is in the +z direction.

.. note::

   While LAMMPS allows specification of a triclinic simulation box in
   either **general** or **restricted** form, internally LAMMPS only
   uses restricted triclinic simulation boxes.  This is for parallel
   efficiency and to formulate partitioning of the simulation box
   across processors, neighbor list building, and inter-processor
   communication of per-atom data with methods similar to those used
   for orthogonal boxes.  This means 3 things.  (1) Input of a general
   triclinic is immediately converted to restricted form.
   (2) If output in general triclinic form is requested (e.g. for atom
   coordinates in a dump file), then conversion from restricted
   triclinic coordinates is done at the time of output.  (3) Most
   importantly, other LAMMPS commands such as the :doc:`region
   <region>` command, that refer to the simulation box geometry,
   operate on restricted triclinic boxes, even if a general triclinic
   box was specified as input.

Note that for 2d simulations a triclinic simulation box is effectively
a parallelogram; see the :doc:'Howto 2d <Howto_2d>` doc page for
details.

The :doc:`boundary <boundary>` command sets the boundary conditions
for the 6 faces of a restricted triclinix box (periodic, non-periodic,
etc), similar to the way the settings apply to the 6 faces of an
orthogonal box.  Note that if a restricted triclinic box is periodic
in the y-dimension and has a non-zero xy tilt factor, then particles
which exit the -y face of the box will re-enter the +y face but will
be displaced in x by the xy tilt factor.  Similarly for z-periodicity,
if the xz and/or yz tilt factors are non-zero, then particles which
exit the -z face of the box will be displaced in x by the xz tilt
factor and in y by the yz tilt factor.

For general and restricted triclinic boxes, their **A**, **B**, **C**
edge vector components can be output via

The :doc:`thermo_style custom <thermo_style>` command also has options
for outputting the parameters that define general and restricted
triclinic simulation boxes.  For general triclinic, this is the
(xlo,ylo,zhi) origin and the 9 components of the **A**, **B**, **C**
edge vectors.  For restricted triclinic, this is (xlo,ylo,zlo),
(xhi,yhi,zhi), and the xy,xz,yz tilt factors.  For both orthogonal and
restricted triclinic boxes, lx/ly/lz refer to the same box sizes,
namely lx = xhi - xlo, etc.

The remainder of this doc page explains mathematical transformations
between different ways of representing general and restrictied
triclinic boxes, which may be useful when creating LAMMPS inputs for
triclinic simulations or interpreting outputs.  How LAMMPS uses tilt
factors for restricted triclinic simulation boxes is also discussed.

----------

Let **A**,\ **B**,\ **C** be the edge vectors of a general triclinic
simulation box.  Assume that **A** x **B** . **C** > 0.  The
equivalent LAMMPS **a**,\ **b**,\ **c** for a restricted triclinic box
are a linear rotation of **A**, **B**, and **C** and can be computed
as follows:

.. math::

  \begin{pmatrix} \mathbf{a}  & \mathbf{b}  & \mathbf{c} \end{pmatrix} = &
  \begin{pmatrix}
    a_x & b_x & c_x \\
    0   & b_y & c_y \\
    0   & 0   & c_z \\
  \end{pmatrix} \\
  a_x = & A \\
  b_x = & \mathbf{B} \cdot \mathbf{\hat{A}} \quad = \quad B \cos{\gamma} \\
  b_y = & |\mathbf{\hat{A}} \times \mathbf{B}| \quad = \quad B \sin{\gamma} \quad =  \quad  \sqrt{B^2 - {b_x}^2} \\
  c_x = & \mathbf{C} \cdot \mathbf{\hat{A}} \quad = \quad C \cos{\beta} \\
  c_y = & \mathbf{C} \cdot \widehat{(\mathbf{A} \times \mathbf{B})} \times \mathbf{\hat{A}} \quad = \quad \frac{\mathbf{B} \cdot \mathbf{C} - b_x c_x}{b_y} \\
  c_z = & |\mathbf{C} \cdot \widehat{(\mathbf{A} \times \mathbf{B})}|\quad = \quad \sqrt{C^2 - {c_x}^2 - {c_y}^2}

where A = \| **A** \| indicates the scalar length of **A**\ . The hat
symbol (\^) indicates the corresponding unit vector. :math:`\beta` and
:math:`\gamma` are angles between the **A**, **B**, **C** vectors
as described below.

If **A** x **B** . **C** < 0 the above equations are not valid for
**c**\ . In this case, it is necessary to first apply an
inversion. This can be achieved by interchanging two of the **A**,
**B**, **C** vectors or by changing the sign of one of them.

For consistency, the same rotation/inversion applied to the triclinic
box edge vectors also typically needs to be applied to atom positions,
velocities, and other vector quantities.  This can be conveniently
achieved by first converting to fractional coordinates in the general
triclinic coordinates and then converting to coordinates in the
resetricted triclinic basis.  The transformation is given by the
following equation:

.. math::

  \mathbf{x} = & \begin{pmatrix} \mathbf{a}  & \mathbf{b}  & \mathbf{c} \end{pmatrix} \cdot \frac{1}{V}
    \begin{pmatrix}
      \mathbf{B \times C}  \\
      \mathbf{C \times A}  \\
      \mathbf{A \times B}
    \end{pmatrix} \cdot \mathbf{X}

where *V* is the volume of the box (same in either basis), **X** is
the fractional vector in the general triclinic basis and **x** is the
resulting vector in the restricted triclinic basis.

----------

General triclinic crystal structures are often defined using three
lattice constants *a*, *b*, and *c*, and three angles :math:`\alpha`,
:math:`\beta`, and :math:`\gamma`. Note that in this nomenclature, the
a, b, and c lattice constants are the scalar lengths of the edge
vectors **a**, **b**, and **c** defined above.  The relationship
between these 6 quantities (a, b, c, :math:`\alpha`, :math:`\beta`,
:math:`\gamma`) and the LAMMPS restricted triclinic box sizes
(lx,ly,lz) = (xhi-xlo,yhi-ylo,zhi-zlo) and tilt factors (xy,xz,yz) is
as follows:

.. math::

  a   = & {\rm lx} \\
  b^2 = &  {\rm ly}^2 +  {\rm xy}^2 \\
  c^2 = &  {\rm lz}^2 +  {\rm xz}^2 +  {\rm yz}^2 \\
  \cos{\alpha} = & \frac{{\rm xy}*{\rm xz} + {\rm ly}*{\rm yz}}{b*c} \\
  \cos{\beta}  = & \frac{\rm xz}{c} \\
  \cos{\gamma} = & \frac{\rm xy}{b} \\

The inverse relationship can be written as follows:

.. math::

  {\rm lx}   = & a \\
  {\rm xy}   = & b \cos{\gamma}  \\
  {\rm xz}   = & c \cos{\beta}\\
  {\rm ly}^2 = & b^2 - {\rm xy}^2 \\
  {\rm yz}   = & \frac{b*c \cos{\alpha} - {\rm xy}*{\rm xz}}{\rm ly} \\
  {\rm lz}^2 = & c^2 - {\rm xz}^2 - {\rm yz}^2 \\

The values of *a*, *b*, *c*, :math:`\alpha` , :math:`\beta`, and
:math:`\gamma` can be printed out or accessed by computes using the
:doc:`thermo_style custom <thermo_style>` keywords *cella*, *cellb*,
*cellc*, *cellalpha*, *cellbeta*, *cellgamma*, respectively.

----------

As discussed on the :doc:`dump <dump>` command doc page, when the BOX
BOUNDS for a snapshot is written to a dump file for a resticted
triclinic box, an orthogonal bounding box which encloses the triclinic
simulation box is output, along with the 3 tilt factors (xy, xz, yz)
of the restricted triclinic box, formatted as follows:

.. parsed-literal::

   ITEM: BOX BOUNDS xy xz yz
   xlo_bound xhi_bound xy
   ylo_bound yhi_bound xz
   zlo_bound zhi_bound yz

This bounding box is convenient for many visualization programs and is
calculated from the 9 restricted triclinic box parameters
(xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz) as follows:

.. parsed-literal::

   xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
   xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
   ylo_bound = ylo + MIN(0.0,yz)
   yhi_bound = yhi + MAX(0.0,yz)
   zlo_bound = zlo
   zhi_bound = zhi

These formulas can be inverted if you need to convert the bounding box
back into the restricted triclinic box parameters, e.g. xlo =
xlo_bound - MIN(0.0,xy,xz,xy+xz).

----------

There is no requirement that a triclinic box be periodic in any
dimension, though as explained above it typically should be in y or z
if you wish enforce a shift in coordinates due to periodic boundary
conditions across the y or z boundaries.

Some commands that work with triclinic boxes, e.g. the :doc:`fix
deform <fix_deform>` and :doc:`fix npt <fix_nh>` commands, require
periodicity or non-shrink-wrap boundary conditions in specific
dimensions.  See the command doc pages for details.

A restricted triclinic box can be defined with all 3 tilt factors =
0.0, so that it is initially orthogonal.  This is necessary if the box
will become non-orthogonal, e.g. due to use of the :doc:`fix npt
<fix_nh>` or :doc:`fix deform <fix_deform>` commands.  Alternatively,
you can use the :doc:`change_box <change_box>` command to convert a
simulation box from orthogonal to restricted triclinic and vice versa.

Highly tilted restricted triclinic simulation boxes can be
computationally inefficient.  This is due to the large volume of
communication needed to acquire ghost atoms around a processor's
irregular-shaped subdomain.  For extreme values of tilt, LAMMPS may
also lose atoms and generate an error.

LAMMPS will issue a warning if you define a restricted triclinic box
with a tilt factor which skews the box more than half the distance of
the parallel box length, which is the first dimension in the tilt
factor (x for xz).

For example, if xlo = 2 and xhi = 12, then the x box length is 10 and
the xy tilt factor should be between -5 and 5 to avoid the warning.
Similarly, both xz and yz should be between -(xhi-xlo)/2 and
+(yhi-ylo)/2.  Note that these are not limitations, since if the
maximum tilt factor is 5 (as in this example), then simulations boxes
and atom configurations with tilt = ..., -15, -5, 5, 15, 25, ... are
geometrically all equivalent.

If the box tilt exceeds this limit during a dynamics run (e.g. due to
the :doc:`fix deform <fix_deform>` command), then by default the box
is "flipped" to an equivalent shape with a tilt factor within the
warning bounds, and the run continues.  See the :doc:`fix deform
<fix_deform>` page for further details.  Box flips that would normally
occur using the :doc:`fix deform <fix_deform>` or :doc:`fix npt
<fix_nh>` commands can be suppressed using the *flip no* option with
either of the commands.

One exception to box flipping is if the first dimension in the tilt
factor (x for xy) is non-periodic.  In that case, the limits on the
tilt factor are not enforced, since flipping the box in that dimension
does not change the atom positions due to non-periodicity.  In this
mode, the system tilts to large angles, the simulation will simply
become inefficient, due to the highly skewed simulation box.

