.. index:: fix deform

fix deform command
==================

fix deform/kk command
=====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID deform N parameter args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* deform = style name of this fix command
* N = perform box deformation every this many timesteps
* one or more parameter/arg pairs may be appended

  .. parsed-literal::

     parameter = *x* or *y* or *z* or *xy* or *xz* or *yz*
       *x*\ , *y*\ , *z* args = style value(s)
         style = *final* or *delta* or *scale* or *vel* or *erate* or *trate* or *volume* or *wiggle* or *variable*
           *final* values = lo hi
             lo hi = box boundaries at end of run (distance units)
           *delta* values = dlo dhi
             dlo dhi = change in box boundaries at end of run (distance units)
           *scale* values = factor
             factor = multiplicative factor for change in box length at end of run
           *vel* value = V
             V = change box length at this velocity (distance/time units),
                 effectively an engineering strain rate
           *erate* value = R
             R = engineering strain rate (1/time units)
           *trate* value = R
             R = true strain rate (1/time units)
           *volume* value = none = adjust this dim to preserve volume of system
           *wiggle* values = A Tp
             A = amplitude of oscillation (distance units)
             Tp = period of oscillation (time units)
           *variable* values = v_name1 v_name2
             v_name1 = variable with name1 for box length change as function of time
             v_name2 = variable with name2 for change rate as function of time
       *xy*\ , *xz*\ , *yz* args = style value
         style = *final* or *delta* or *vel* or *erate* or *trate* or *wiggle*
           *final* value = tilt
             tilt = tilt factor at end of run (distance units)
           *delta* value = dtilt
             dtilt = change in tilt factor at end of run (distance units)
           *vel* value = V
             V = change tilt factor at this velocity (distance/time units),
                 effectively an engineering shear strain rate
           *erate* value = R
             R = engineering shear strain rate (1/time units)
           *trate* value = R
             R = true shear strain rate (1/time units)
           *wiggle* values = A Tp
             A = amplitude of oscillation (distance units)
             Tp = period of oscillation (time units)
           *variable* values = v_name1 v_name2
             v_name1 = variable with name1 for tilt change as function of time
             v_name2 = variable with name2 for change rate as function of time

* zero or more keyword/value pairs may be appended
* keyword = *remap* or *flip* or *units*

  .. parsed-literal::

       *remap* value = *x* or *v* or *none*
         x = remap coords of atoms in group into deforming box
         v = remap velocities of all atoms when they cross periodic boundaries
         none = no remapping of x or v
       *flip* value = *yes* or *no*
         allow or disallow box flips when it becomes highly skewed
       *units* value = *lattice* or *box*
         lattice = distances are defined in lattice units
         box = distances are defined in simulation box units

Examples
""""""""

.. parsed-literal::

   fix 1 all deform 1 x final 0.0 9.0 z final 0.0 5.0 units box
   fix 1 all deform 1 x trate 0.1 y volume z volume
   fix 1 all deform 1 xy erate 0.001 remap v
   fix 1 all deform 10 y delta -0.5 0.5 xz vel 1.0

Description
"""""""""""

Change the volume and/or shape of the simulation box during a dynamics
run.  Orthogonal simulation boxes have 3 adjustable parameters
(x,y,z).  Triclinic (non-orthogonal) simulation boxes have 6
adjustable parameters (x,y,z,xy,xz,yz).  Any or all of them can be
adjusted independently and simultaneously by this command.

This fix can be used to perform non-equilibrium MD (NEMD) simulations
of a continuously strained system.  See the :doc:`fix nvt/sllod <fix_nvt_sllod>` and :doc:`compute temp/deform <compute_temp_deform>` commands for more details.  Note
that simulation of a continuously extended system (extensional flow)
can be modeled using the :ref:`USER-UEF package <PKG-USER-UEF>` and its :doc:`fix commands <fix_nh_uef>`.

For the *x*\ , *y*\ , *z* parameters, the associated dimension cannot be
shrink-wrapped.  For the *xy*\ , *yz*\ , *xz* parameters, the associated
2nd dimension cannot be shrink-wrapped.  Dimensions not varied by this
command can be periodic or non-periodic.  Dimensions corresponding to
unspecified parameters can also be controlled by a :doc:`fix npt <fix_nh>` or :doc:`fix nph <fix_nh>` command.

The size and shape of the simulation box at the beginning of the
simulation run were either specified by the
:doc:`create_box <create_box>` or :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` command used to setup the simulation
initially if it is the first run, or they are the values from the end
of the previous run.  The :doc:`create_box <create_box>`, :doc:`read data <read_data>`, and :doc:`read_restart <read_restart>` commands
specify whether the simulation box is orthogonal or non-orthogonal
(triclinic) and explain the meaning of the xy,xz,yz tilt factors.  If
fix deform changes the xy,xz,yz tilt factors, then the simulation box
must be triclinic, even if its initial tilt factors are 0.0.

As described below, the desired simulation box size and shape at the
end of the run are determined by the parameters of the fix deform
command.  Every Nth timestep during the run, the simulation box is
expanded, contracted, or tilted to ramped values between the initial
and final values.

----------

For the *x*\ , *y*\ , and *z* parameters, this is the meaning of their
styles and values.

The *final*\ , *delta*\ , *scale*\ , *vel*\ , and *erate* styles all change
the specified dimension of the box via "constant displacement" which
is effectively a "constant engineering strain rate".  This means the
box dimension changes linearly with time from its initial to final
value.

For style *final*\ , the final lo and hi box boundaries of a dimension
are specified.  The values can be in lattice or box distance units.
See the discussion of the units keyword below.

For style *delta*\ , plus or minus changes in the lo/hi box boundaries
of a dimension are specified.  The values can be in lattice or box
distance units.  See the discussion of the units keyword below.

For style *scale*\ , a multiplicative factor to apply to the box length
of a dimension is specified.  For example, if the initial box length
is 10, and the factor is 1.1, then the final box length will be 11.  A
factor less than 1.0 means compression.

For style *vel*\ , a velocity at which the box length changes is
specified in units of distance/time.  This is effectively a "constant
engineering strain rate", where rate = V/L0 and L0 is the initial box
length.  The distance can be in lattice or box distance units.  See
the discussion of the units keyword below.  For example, if the
initial box length is 100 Angstroms, and V is 10 Angstroms/psec, then
after 10 psec, the box length will have doubled.  After 20 psec, it
will have tripled.

The *erate* style changes a dimension of the box at a "constant
engineering strain rate".  The units of the specified strain rate are
1/time.  See the :doc:`units <units>` command for the time units
associated with different choices of simulation units,
e.g. picoseconds for "metal" units).  Tensile strain is unitless and
is defined as delta/L0, where L0 is the original box length and delta
is the change relative to the original length.  The box length L as a
function of time will change as

.. parsed-literal::

   L(t) = L0 (1 + erate\*dt)

where dt is the elapsed time (in time units).  Thus if *erate* R is
specified as 0.1 and time units are picoseconds, this means the box
length will increase by 10% of its original length every picosecond.
I.e. strain after 1 psec = 0.1, strain after 2 psec = 0.2, etc.  R =
-0.01 means the box length will shrink by 1% of its original length
every picosecond.  Note that for an "engineering" rate the change is
based on the original box length, so running with R = 1 for 10
picoseconds expands the box length by a factor of 11 (strain of 10),
which is different that what the *trate* style would induce.

The *trate* style changes a dimension of the box at a "constant true
strain rate".  Note that this is not an "engineering strain rate", as
the other styles are.  Rather, for a "true" rate, the rate of change
is constant, which means the box dimension changes non-linearly with
time from its initial to final value.  The units of the specified
strain rate are 1/time.  See the :doc:`units <units>` command for the
time units associated with different choices of simulation units,
e.g. picoseconds for "metal" units).  Tensile strain is unitless and
is defined as delta/L0, where L0 is the original box length and delta
is the change relative to the original length.

The box length L as a function of time will change as

.. parsed-literal::

   L(t) = L0 exp(trate\*dt)

where dt is the elapsed time (in time units).  Thus if *trate* R is
specified as ln(1.1) and time units are picoseconds, this means the
box length will increase by 10% of its current (not original) length
every picosecond.  I.e. strain after 1 psec = 0.1, strain after 2 psec
= 0.21, etc.  R = ln(2) or ln(3) means the box length will double or
triple every picosecond.  R = ln(0.99) means the box length will
shrink by 1% of its current length every picosecond.  Note that for a
"true" rate the change is continuous and based on the current length,
so running with R = ln(2) for 10 picoseconds does not expand the box
length by a factor of 11 as it would with *erate*\ , but by a factor of
1024 since the box length will double every picosecond.

Note that to change the volume (or cross-sectional area) of the
simulation box at a constant rate, you can change multiple dimensions
via *erate* or *trate*\ .  E.g. to double the box volume in a picosecond
picosecond, you could set "x erate M", "y erate M", "z erate M", with
M = pow(2,1/3) - 1 = 0.26, since if each box dimension grows by 26%,
the box volume doubles.  Or you could set "x trate M", "y trate M", "z
trate M", with M = ln(1.26) = 0.231, and the box volume would double
every picosecond.

The *volume* style changes the specified dimension in such a way that
the box volume remains constant while other box dimensions are changed
explicitly via the styles discussed above.  For example, "x scale 1.1
y scale 1.1 z volume" will shrink the z box length as the x,y box
lengths increase, to keep the volume constant (product of x,y,z
lengths).  If "x scale 1.1 z volume" is specified and parameter *y* is
unspecified, then the z box length will shrink as x increases to keep
the product of x,z lengths constant.  If "x scale 1.1 y volume z
volume" is specified, then both the y,z box lengths will shrink as x
increases to keep the volume constant (product of x,y,z lengths).  In
this case, the y,z box lengths shrink so as to keep their relative
aspect ratio constant.

For solids or liquids, note that when one dimension of the box is
expanded via fix deform (i.e. tensile strain), it may be physically
undesirable to hold the other 2 box lengths constant (unspecified by
fix deform) since that implies a density change.  Using the *volume*
style for those 2 dimensions to keep the box volume constant may make
more physical sense, but may also not be correct for materials and
potentials whose Poisson ratio is not 0.5.  An alternative is to use
:doc:`fix npt aniso <fix_nh>` with zero applied pressure on those 2
dimensions, so that they respond to the tensile strain dynamically.

The *wiggle* style oscillates the specified box length dimension
sinusoidally with the specified amplitude and period.  I.e. the box
length L as a function of time is given by

.. parsed-literal::

   L(t) = L0 + A sin(2\*pi t/Tp)

where L0 is its initial length.  If the amplitude A is a positive
number the box initially expands, then contracts, etc.  If A is
negative then the box initially contracts, then expands, etc.  The
amplitude can be in lattice or box distance units.  See the discussion
of the units keyword below.

The *variable* style changes the specified box length dimension by
evaluating a variable, which presumably is a function of time.  The
variable with *name1* must be an :doc:`equal-style variable <variable>`
and should calculate a change in box length in units of distance.
Note that this distance is in box units, not lattice units; see the
discussion of the *units* keyword below.  The formula associated with
variable *name1* can reference the current timestep.  Note that it
should return the "change" in box length, not the absolute box length.
This means it should evaluate to 0.0 when invoked on the initial
timestep of the run following the definition of fix deform.  It should
evaluate to a value > 0.0 to dilate the box at future times, or a
value < 0.0 to compress the box.

The variable *name2* must also be an :doc:`equal-style variable <variable>` and should calculate the rate of box length
change, in units of distance/time, i.e. the time-derivative of the
*name1* variable.  This quantity is used internally by LAMMPS to reset
atom velocities when they cross periodic boundaries.  It is computed
internally for the other styles, but you must provide it when using an
arbitrary variable.

Here is an example of using the *variable* style to perform the same
box deformation as the *wiggle* style formula listed above, where we
assume that the current timestep = 0.

.. parsed-literal::

   variable A equal 5.0
   variable Tp equal 10.0
   variable displace equal "v_A \* sin(2\*PI \* step\*dt/v_Tp)"
   variable rate equal "2\*PI\*v_A/v_Tp \* cos(2\*PI \* step\*dt/v_Tp)"
   fix 2 all deform 1 x variable v_displace v_rate remap v

For the *scale*\ , *vel*\ , *erate*\ , *trate*\ , *volume*\ , *wiggle*\ , and
*variable* styles, the box length is expanded or compressed around its
mid point.

----------

For the *xy*\ , *xz*\ , and *yz* parameters, this is the meaning of their
styles and values.  Note that changing the tilt factors of a triclinic
box does not change its volume.

The *final*\ , *delta*\ , *vel*\ , and *erate* styles all change the shear
strain at a "constant engineering shear strain rate".  This means the
tilt factor changes linearly with time from its initial to final
value.

For style *final*\ , the final tilt factor is specified.  The value
can be in lattice or box distance units.  See the discussion of the
units keyword below.

For style *delta*\ , a plus or minus change in the tilt factor is
specified.  The value can be in lattice or box distance units.  See
the discussion of the units keyword below.

For style *vel*\ , a velocity at which the tilt factor changes is
specified in units of distance/time.  This is effectively an
"engineering shear strain rate", where rate = V/L0 and L0 is the
initial box length perpendicular to the direction of shear.  The
distance can be in lattice or box distance units.  See the discussion
of the units keyword below.  For example, if the initial tilt factor
is 5 Angstroms, and the V is 10 Angstroms/psec, then after 1 psec, the
tilt factor will be 15 Angstroms.  After 2 psec, it will be 25
Angstroms.

The *erate* style changes a tilt factor at a "constant engineering
shear strain rate".  The units of the specified shear strain rate are
1/time.  See the :doc:`units <units>` command for the time units
associated with different choices of simulation units,
e.g. picoseconds for "metal" units).  Shear strain is unitless and is
defined as offset/length, where length is the box length perpendicular
to the shear direction (e.g. y box length for xy deformation) and
offset is the displacement distance in the shear direction (e.g. x
direction for xy deformation) from the unstrained orientation.

The tilt factor T as a function of time will change as

.. parsed-literal::

   T(t) = T0 + L0\*erate\*dt

where T0 is the initial tilt factor, L0 is the original length of the
box perpendicular to the shear direction (e.g. y box length for xy
deformation), and dt is the elapsed time (in time units).  Thus if
*erate* R is specified as 0.1 and time units are picoseconds, this
means the shear strain will increase by 0.1 every picosecond.  I.e. if
the xy shear strain was initially 0.0, then strain after 1 psec = 0.1,
strain after 2 psec = 0.2, etc.  Thus the tilt factor would be 0.0 at
time 0, 0.1\*ybox at 1 psec, 0.2\*ybox at 2 psec, etc, where ybox is the
original y box length.  R = 1 or 2 means the tilt factor will increase
by 1 or 2 every picosecond.  R = -0.01 means a decrease in shear
strain by 0.01 every picosecond.

The *trate* style changes a tilt factor at a "constant true shear
strain rate".  Note that this is not an "engineering shear strain
rate", as the other styles are.  Rather, for a "true" rate, the rate
of change is constant, which means the tilt factor changes
non-linearly with time from its initial to final value.  The units of
the specified shear strain rate are 1/time.  See the
:doc:`units <units>` command for the time units associated with
different choices of simulation units, e.g. picoseconds for "metal"
units).  Shear strain is unitless and is defined as offset/length,
where length is the box length perpendicular to the shear direction
(e.g. y box length for xy deformation) and offset is the displacement
distance in the shear direction (e.g. x direction for xy deformation)
from the unstrained orientation.

The tilt factor T as a function of time will change as

.. parsed-literal::

   T(t) = T0 exp(trate\*dt)

where T0 is the initial tilt factor and dt is the elapsed time (in
time units).  Thus if *trate* R is specified as ln(1.1) and time units
are picoseconds, this means the shear strain or tilt factor will
increase by 10% every picosecond.  I.e. if the xy shear strain was
initially 0.1, then strain after 1 psec = 0.11, strain after 2 psec =
0.121, etc.  R = ln(2) or ln(3) means the tilt factor will double or
triple every picosecond.  R = ln(0.99) means the tilt factor will
shrink by 1% every picosecond.  Note that the change is continuous, so
running with R = ln(2) for 10 picoseconds does not change the tilt
factor by a factor of 10, but by a factor of 1024 since it doubles
every picosecond.  Note that the initial tilt factor must be non-zero
to use the *trate* option.

Note that shear strain is defined as the tilt factor divided by the
perpendicular box length.  The *erate* and *trate* styles control the
tilt factor, but assume the perpendicular box length remains constant.
If this is not the case (e.g. it changes due to another fix deform
parameter), then this effect on the shear strain is ignored.

The *wiggle* style oscillates the specified tilt factor sinusoidally
with the specified amplitude and period.  I.e. the tilt factor T as a
function of time is given by

.. parsed-literal::

   T(t) = T0 + A sin(2\*pi t/Tp)

where T0 is its initial value.  If the amplitude A is a positive
number the tilt factor initially becomes more positive, then more
negative, etc.  If A is negative then the tilt factor initially
becomes more negative, then more positive, etc.  The amplitude can be
in lattice or box distance units.  See the discussion of the units
keyword below.

The *variable* style changes the specified tilt factor by evaluating a
variable, which presumably is a function of time.  The variable with
*name1* must be an :doc:`equal-style variable <variable>` and should
calculate a change in tilt in units of distance.  Note that this
distance is in box units, not lattice units; see the discussion of the
*units* keyword below.  The formula associated with variable *name1*
can reference the current timestep.  Note that it should return the
"change" in tilt factor, not the absolute tilt factor.  This means it
should evaluate to 0.0 when invoked on the initial timestep of the run
following the definition of fix deform.

The variable *name2* must also be an :doc:`equal-style variable <variable>` and should calculate the rate of tilt change,
in units of distance/time, i.e. the time-derivative of the *name1*
variable.  This quantity is used internally by LAMMPS to reset atom
velocities when they cross periodic boundaries.  It is computed
internally for the other styles, but you must provide it when using an
arbitrary variable.

Here is an example of using the *variable* style to perform the same
box deformation as the *wiggle* style formula listed above, where we
assume that the current timestep = 0.

.. parsed-literal::

   variable A equal 5.0
   variable Tp equal 10.0
   variable displace equal "v_A \* sin(2\*PI \* step\*dt/v_Tp)"
   variable rate equal "2\*PI\*v_A/v_Tp \* cos(2\*PI \* step\*dt/v_Tp)"
   fix 2 all deform 1 xy variable v_displace v_rate remap v

----------

All of the tilt styles change the xy, xz, yz tilt factors during a
simulation.  In LAMMPS, tilt factors (xy,xz,yz) for triclinic boxes
are normally bounded by half the distance of the parallel box length.
See the discussion of the *flip* keyword below, to allow this bound to
be exceeded, if desired.

For example, if xlo = 2 and xhi = 12, then the x box length is 10 and
the xy tilt factor must be between -5 and 5.  Similarly, both xz and
yz must be between -(xhi-xlo)/2 and +(yhi-ylo)/2.  Note that this is
not a limitation, since if the maximum tilt factor is 5 (as in this
example), then configurations with tilt = ..., -15, -5, 5, 15, 25,
... are all equivalent.

To obey this constraint and allow for large shear deformations to be
applied via the *xy*\ , *xz*\ , or *yz* parameters, the following
algorithm is used.  If *prd* is the associated parallel box length (10
in the example above), then if the tilt factor exceeds the accepted
range of -5 to 5 during the simulation, then the box is flipped to the
other limit (an equivalent box) and the simulation continues.  Thus
for this example, if the initial xy tilt factor was 0.0 and "xy final
100.0" was specified, then during the simulation the xy tilt factor
would increase from 0.0 to 5.0, the box would be flipped so that the
tilt factor becomes -5.0, the tilt factor would increase from -5.0 to
5.0, the box would be flipped again, etc.  The flip occurs 10 times
and the final tilt factor at the end of the simulation would be 0.0.
During each flip event, atoms are remapped into the new box in the
appropriate manner.

The one exception to this rule is if the 1st dimension in the tilt
factor (x for xy) is non-periodic.  In that case, the limits on the
tilt factor are not enforced, since flipping the box in that dimension
does not change the atom positions due to non-periodicity.  In this
mode, if you tilt the system to extreme angles, the simulation will
simply become inefficient due to the highly skewed simulation box.

----------

Each time the box size or shape is changed, the *remap* keyword
determines whether atom positions are remapped to the new box.  If
*remap* is set to *x* (the default), atoms in the fix group are
remapped; otherwise they are not.  Note that their velocities are not
changed, just their positions are altered.  If *remap* is set to *v*\ ,
then any atom in the fix group that crosses a periodic boundary will
have a delta added to its velocity equal to the difference in
velocities between the lo and hi boundaries.  Note that this velocity
difference can include tilt components, e.g. a delta in the x velocity
when an atom crosses the y periodic boundary.  If *remap* is set to
*none*\ , then neither of these remappings take place.

Conceptually, setting *remap* to *x* forces the atoms to deform via an
affine transformation that exactly matches the box deformation.  This
setting is typically appropriate for solids.  Note that though the
atoms are effectively "moving" with the box over time, it is not due
to their having a velocity that tracks the box change, but only due to
the remapping.  By contrast, setting *remap* to *v* is typically
appropriate for fluids, where you want the atoms to respond to the
change in box size/shape on their own and acquire a velocity that
matches the box change, so that their motion will naturally track the
box without explicit remapping of their coordinates.

.. note::

   When non-equilibrium MD (NEMD) simulations are performed using
   this fix, the option "remap v" should normally be used.  This is
   because :doc:`fix nvt/sllod <fix_nvt_sllod>` adjusts the atom positions
   and velocities to induce a velocity profile that matches the changing
   box size/shape.  Thus atom coordinates should NOT be remapped by fix
   deform, but velocities SHOULD be when atoms cross periodic boundaries,
   since that is consistent with maintaining the velocity profile already
   created by fix nvt/sllod.  LAMMPS will warn you if the *remap* setting
   is not consistent with fix nvt/sllod.

.. note::

   For non-equilibrium MD (NEMD) simulations using "remap v" it is
   usually desirable that the fluid (or flowing material, e.g. granular
   particles) stream with a velocity profile consistent with the
   deforming box.  As mentioned above, using a thermostat such as :doc:`fix nvt/sllod <fix_nvt_sllod>` or :doc:`fix lavgevin <fix_langevin>`
   (with a bias provided by :doc:`compute temp/deform <compute_temp_deform>`), will typically accomplish
   that.  If you do not use a thermostat, then there is no driving force
   pushing the atoms to flow in a manner consistent with the deforming
   box.  E.g. for a shearing system the box deformation velocity may vary
   from 0 at the bottom to 10 at the top of the box.  But the stream
   velocity profile of the atoms may vary from -5 at the bottom to +5 at
   the top.  You can monitor these effects using the :doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`compute temp/deform <compute_temp_deform>`, and :doc:`compute temp/profile <compute_temp_profile>` commands.  One way to induce
   atoms to stream consistent with the box deformation is to give them an
   initial velocity profile, via the :doc:`velocity ramp <velocity>`
   command, that matches the box deformation rate.  This also typically
   helps the system come to equilibrium more quickly, even if a
   thermostat is used.

.. note::

   If a :doc:`fix rigid <fix_rigid>` is defined for rigid bodies, and
   *remap* is set to *x*\ , then the center-of-mass coordinates of rigid
   bodies will be remapped to the changing simulation box.  This will be
   done regardless of whether atoms in the rigid bodies are in the fix
   deform group or not.  The velocity of the centers of mass are not
   remapped even if *remap* is set to *v*\ , since :doc:`fix nvt/sllod <fix_nvt_sllod>` does not currently do anything special
   for rigid particles.  If you wish to perform a NEMD simulation of
   rigid particles, you can either thermostat them independently or
   include a background fluid and thermostat the fluid via :doc:`fix nvt/sllod <fix_nvt_sllod>`.

The *flip* keyword allows the tilt factors for a triclinic box to
exceed half the distance of the parallel box length, as discussed
above.  If the *flip* value is set to *yes*\ , the bound is enforced by
flipping the box when it is exceeded.  If the *flip* value is set to
*no*\ , the tilt will continue to change without flipping.  Note that if
you apply large deformations, this means the box shape can tilt
dramatically LAMMPS will run less efficiently, due to the large volume
of communication needed to acquire ghost atoms around a processor's
irregular-shaped sub-domain.  For extreme values of tilt, LAMMPS may
also lose atoms and generate an error.

The *units* keyword determines the meaning of the distance units used
to define various arguments.  A *box* value selects standard distance
units as defined by the :doc:`units <units>` command, e.g. Angstroms for
units = real or metal.  A *lattice* value means the distance units are
in lattice spacings.  The :doc:`lattice <lattice>` command must have
been previously used to define the lattice spacing.  Note that the
units choice also affects the *vel* style parameters since it is
defined in terms of distance/time.  Also note that the units keyword
does not affect the *variable* style.  You should use the *xlat*\ ,
*ylat*\ , *zlat* keywords of the :doc:`thermo_style <thermo_style>`
command if you want to include lattice spacings in a variable formula.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix will restore the initial box settings from :doc:`binary restart files <restart>`, which allows the fix to be properly continue
deformation, when using the start/stop options of the :doc:`run <run>`
command.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.

This fix can perform deformation over multiple runs, using the *start*
and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

You cannot apply x, y, or z deformations to a dimension that is
shrink-wrapped via the :doc:`boundary <boundary>` command.

You cannot apply xy, yz, or xz deformations to a 2nd dimension (y in
xy) that is shrink-wrapped via the :doc:`boundary <boundary>` command.

Related commands
""""""""""""""""

:doc:`change_box <change_box>`

Default
"""""""

The option defaults are remap = x, flip = yes, and units = lattice.
