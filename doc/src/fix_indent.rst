.. index:: fix indent

fix indent command
==================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID indent K keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* indent = style name of this fix command
* K = force constant for indenter surface (force/distance\^2 units)
* one or more keyword/value pairs may be appended
* keyword = *sphere* or *cylinder* or *plane* or *side* or *units*

  .. parsed-literal::

       *sphere* args = x y z R
         x,y,z = initial position of center of indenter (distance units)
         R = sphere radius of indenter (distance units)
         any of x,y,z,R can be a variable (see below)
       *cylinder* args = dim c1 c2 R
         dim = *x* or *y* or *z* = axis of cylinder
         c1,c2 = coords of cylinder axis in other 2 dimensions (distance units)
         R = cylinder radius of indenter (distance units)
         any of c1,c2,R can be a variable (see below)
       *plane* args = dim pos side
         dim = *x* or *y* or *z* = plane perpendicular to this dimension
         pos = position of plane in dimension x, y, or z (distance units)
         pos can be a variable (see below)
         side = *lo* or *hi*
       *side* value = *in* or *out*
         *in* = the indenter acts on particles inside the sphere or cylinder
         *out* = the indenter acts on particles outside the sphere or cylinder
       *units* value = *lattice* or *box*
         lattice = the geometry is defined in lattice units
         box = the geometry is defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all indent 10.0 sphere 0.0 0.0 15.0 3.0
   fix 1 all indent 10.0 sphere v_x v_y 0.0 v_radius side in
   fix 2 flow indent 10.0 cylinder z 0.0 0.0 10.0 units box

Description
"""""""""""

Insert an indenter within a simulation box.  The indenter repels all
atoms in the group that touch it, so it can be used to push into a
material or as an obstacle in a flow.  Or it can be used as a
constraining wall around a simulation; see the discussion of the
*side* keyword below.

The indenter can either be spherical or cylindrical or planar.  You
must set one of those 3 keywords.

A spherical indenter exerts a force of magnitude

.. math::

   F(r) = - K \left( r - R \right)^2

on each atom where *K* is the specified force constant, *r* is the
distance from the atom to the center of the indenter, and *R* is the
radius of the indenter.  The force is repulsive and F(r) = 0 for *r* >
*R*\ .

A cylindrical indenter exerts the same force, except that *r* is the
distance from the atom to the center axis of the cylinder.  The
cylinder extends infinitely along its axis.

Spherical and cylindrical indenters account for periodic boundaries in
two ways.  First, the center point of a spherical indenter (x,y,z) or
axis of a cylindrical indenter (c1,c2) is remapped back into the
simulation box, if the box is periodic in a particular dimension.
This occurs every timestep if the indenter geometry is specified with
a variable (see below), e.g. it is moving over time.  Second, the
calculation of distance to the indenter center or axis accounts for
periodic boundaries.  Both of these mean that an indenter can
effectively move through and straddle one or more periodic boundaries.

A planar indenter is really an axis-aligned infinite-extent wall
exerting the same force on atoms in the system, where *R* is the
position of the plane and *r-R* is the distance from the plane.  If
the *side* parameter of the plane is specified as *lo* then it will
indent from the lo end of the simulation box, meaning that atoms with
a coordinate less than the plane's current position will be pushed
towards the hi end of the box and atoms with a coordinate higher than
the plane's current position will feel no force.  Vice versa if *side*
is specified as *hi*\ .

Any of the 4 quantities defining a spherical indenter's geometry can
be specified as an equal-style :doc:`variable <variable>`, namely *x*\ ,
*y*\ , *z*\ , or *R*\ .  Similarly, for a cylindrical indenter, any of *c1*\ ,
*c2*\ , or *R*\ , can be a variable.  For a planar indenter, *pos* can be
a variable.  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to define the
indenter geometry.

Note that equal-style variables can specify formulas with various
mathematical functions, and include :doc:`thermo_style <thermo_style>`
command keywords for the simulation box parameters and timestep and
elapsed time.  Thus it is easy to specify indenter properties that
change as a function of time or span consecutive runs in a continuous
fashion.  For the latter, see the *start* and *stop* keywords of the
:doc:`run <run>` command and the *elaplong* keyword of :doc:`thermo_style custom <thermo_style>` for details.

For example, if a spherical indenter's x-position is specified as v\_x,
then this variable definition will keep it's center at a relative
position in the simulation box, 1/4 of the way from the left edge to
the right edge, even if the box size changes:

.. code-block:: LAMMPS

   variable x equal "xlo + 0.25*lx"

Similarly, either of these variable definitions will move the indenter
from an initial position at 2.5 at a constant velocity of 5:

.. code-block:: LAMMPS

   variable x equal "2.5 + 5*elaplong*dt"
   variable x equal vdisplace(2.5,5)

If a spherical indenter's radius is specified as v\_r, then these
variable definitions will grow the size of the indenter at a specified
rate.

.. code-block:: LAMMPS

   variable r0 equal 0.0
   variable rate equal 1.0
   variable r equal "v_r0 + step*dt*v_rate"

If the *side* keyword is specified as *out*\ , which is the default,
then particles outside the indenter are pushed away from its outer
surface, as described above.  This only applies to spherical or
cylindrical indenters.  If the *side* keyword is specified as *in*\ ,
the action of the indenter is reversed.  Particles inside the indenter
are pushed away from its inner surface.  In other words, the indenter
is now a containing wall that traps the particles inside it.  If the
radius shrinks over time, it will squeeze the particles.

The *units* keyword determines the meaning of the distance units used
to define the indenter geometry.  A *box* value selects standard
distance units as defined by the :doc:`units <units>` command,
e.g. Angstroms for units = real or metal.  A *lattice* value means the
distance units are in lattice spacings.  The :doc:`lattice <lattice>`
command must have been previously used to define the lattice spacing.
The (x,y,z) coords of the indenter position are scaled by the x,y,z
lattice spacings respectively.  The radius of a spherical or
cylindrical indenter is scaled by the x lattice spacing.

Note that the units keyword only affects indenter geometry parameters
specified directly with numbers, not those specified as variables.  In
the latter case, you should use the *xlat*\ , *ylat*\ , *zlat* keywords of
the :doc:`thermo_style <thermo_style>` command if you want to include
lattice spacings in a variable formula.

The force constant *K* is not affected by the *units* keyword.  It is
always in force/distance\^2 units where force and distance are defined
by the :doc:`units <units>` command.  If you wish K to be scaled by the
lattice spacing, you can define K with a variable whose formula
contains *xlat*\ , *ylat*\ , *zlat* keywords of the
:doc:`thermo_style <thermo_style>` command, e.g.

.. code-block:: LAMMPS

   variable k equal 100.0/xlat/xlat
   fix 1 all indent $k sphere ...

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy of interaction between atoms and the indenter to
the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.  The energy of each particle interacting
with the indenter is K/3 (r - R)\^3.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar energy and a global 3-vector of
forces (on the indenter), which can be accessed by various :doc:`output commands <Howto_output>`.  The scalar and vector values calculated
by this fix are "extensive".

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  Note that if you
define the indenter geometry with a variable using a time-dependent
formula, LAMMPS uses the iteration count in the minimizer as the
timestep.  But it is almost certainly a bad idea to have the indenter
change its position or size during a minimization.  LAMMPS does not
check if you have done this.

.. note::

   If you want the atom/indenter interaction energy to be included
   in the total potential energy of the system (the quantity being
   minimized), you must enable the :doc:`fix_modify <fix_modify>` *energy*
   option for this fix.

Restrictions
""""""""""""
 none

**Related commands:** none

Default
"""""""

The option defaults are side = out and units = lattice.
