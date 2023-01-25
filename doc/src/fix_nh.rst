.. index:: fix nvt
.. index:: fix nvt/gpu
.. index:: fix nvt/intel
.. index:: fix nvt/kk
.. index:: fix nvt/omp
.. index:: fix npt
.. index:: fix npt/gpu
.. index:: fix npt/intel
.. index:: fix npt/kk
.. index:: fix npt/omp
.. index:: fix nph
.. index:: fix nph/kk
.. index:: fix nph/omp

fix nvt command
===============

Accelerator Variants: *nvt/gpu*, *nvt/intel*, *nvt/kk*, *nvt/omp*

fix npt command
===============

Accelerator Variants: *npt/gpu*, *npt/intel*, *npt/kk*, *npt/omp*

fix nph command
===============

Accelerator Variants: *nph/kk*, *nph/omp*

Syntax
""""""

.. parsed-literal::

   fix ID group-ID style_name keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *nvt* or *npt* or *nph*
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *temp* or *iso* or *aniso* or *tri* or *x* or *y* or *z* or *xy* or *yz* or *xz* or *couple* or *tchain* or *pchain* or *mtk* or *tloop* or *ploop* or *nreset* or *drag* or *ptemp* or *dilate* or *scalexy* or *scaleyz* or *scalexz* or *flip* or *fixedpoint* or *update*
       *temp* values = Tstart Tstop Tdamp
         Tstart,Tstop = external temperature at start/end of run
         Tdamp = temperature damping parameter (time units)
       *iso* or *aniso* or *tri* values = Pstart Pstop Pdamp
         Pstart,Pstop = scalar external pressure at start/end of run (pressure units)
         Pdamp = pressure damping parameter (time units)
       *x* or *y* or *z* or *xy* or *yz* or *xz* values = Pstart Pstop Pdamp
         Pstart,Pstop = external stress tensor component at start/end of run (pressure units)
         Pdamp = stress damping parameter (time units)
       *couple* = *none* or *xyz* or *xy* or *yz* or *xz*
       *tchain* value = N
         N = length of thermostat chain (1 = single thermostat)
       *pchain* value = N
         N length of thermostat chain on barostat (0 = no thermostat)
       *mtk* value = *yes* or *no* = add in MTK adjustment term or not
       *tloop* value = M
         M = number of sub-cycles to perform on thermostat
       *ploop* value = M
         M = number of sub-cycles to perform on barostat thermostat
       *nreset* value = reset reference cell every this many timesteps
       *drag* value = Df
         Df = drag factor added to barostat/thermostat (0.0 = no drag)
       *ptemp* value = Ttarget
         Ttarget = target temperature for barostat
       *dilate* value = dilate-group-ID
         dilate-group-ID = only dilate atoms in this group due to barostat volume changes
       *scalexy* value = *yes* or *no* = scale xy with ly
       *scaleyz* value = *yes* or *no* = scale yz with lz
       *scalexz* value = *yes* or *no* = scale xz with lz
       *flip* value = *yes* or *no* = allow or disallow box flips when it becomes highly skewed
       *fixedpoint* values = x y z
         x,y,z = perform barostat dilation/contraction around this point (distance units)
       *update* value = *dipole* or *dipole/dlm*
         dipole = update dipole orientation (only for sphere variants)
         dipole/dlm = use DLM integrator to update dipole orientation (only for sphere variants)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt temp 300.0 300.0 100.0
   fix 1 water npt temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0
   fix 2 jello npt temp 300.0 300.0 100.0 tri 5.0 5.0 1000.0
   fix 2 ice nph x 1.0 1.0 0.5 y 2.0 2.0 0.5 z 3.0 3.0 0.5 yz 0.1 0.1 0.5 xz 0.2 0.2 0.5 xy 0.3 0.3 0.5 nreset 1000

Description
"""""""""""

These commands perform time integration on Nose-Hoover style
non-Hamiltonian equations of motion which are designed to generate
positions and velocities sampled from the canonical (nvt),
isothermal-isobaric (npt), and isenthalpic (nph) ensembles.  This
updates the position and velocity for atoms in the group each
timestep.

The thermostatting and barostatting is achieved by adding some dynamic
variables which are coupled to the particle velocities
(thermostatting) and simulation domain dimensions (barostatting).  In
addition to basic thermostatting and barostatting, these fixes can
also create a chain of thermostats coupled to the particle thermostat,
and another chain of thermostats coupled to the barostat
variables. The barostat can be coupled to the overall box volume, or
to individual dimensions, including the *xy*, *xz* and *yz* tilt
dimensions. The external pressure of the barostat can be specified as
either a scalar pressure (isobaric ensemble) or as components of a
symmetric stress tensor (constant stress ensemble).  When used
correctly, the time-averaged temperature and stress tensor of the
particles will match the target values specified by Tstart/Tstop and
Pstart/Pstop.

The equations of motion used are those of Shinoda et al in
:ref:`(Shinoda) <nh-Shinoda>`, which combine the hydrostatic equations of
Martyna, Tobias and Klein in :ref:`(Martyna) <nh-Martyna>` with the strain
energy proposed by Parrinello and Rahman in
:ref:`(Parrinello) <nh-Parrinello>`.  The time integration schemes closely
follow the time-reversible measure-preserving Verlet and rRESPA
integrators derived by Tuckerman et al in :ref:`(Tuckerman) <nh-Tuckerman>`.

----------

The thermostat parameters for fix styles *nvt* and *npt* are specified
using the *temp* keyword.  Other thermostat-related keywords are
*tchain*, *tloop* and *drag*, which are discussed below.

The thermostat is applied to only the translational degrees of freedom
for the particles.  The translational degrees of freedom can also have
a bias velocity removed before thermostatting takes place; see the
description below.  The desired temperature at each timestep is a
ramped value during the run from *Tstart* to *Tstop*\ .  The *Tdamp*
parameter is specified in time units and determines how rapidly the
temperature is relaxed.  For example, a value of 10.0 means to relax
the temperature in a timespan of (roughly) 10 time units (e.g. :math:`\tau`
or fs or ps - see the :doc:`units <units>` command).  The atoms in the
fix group are the only ones whose velocities and positions are updated
by the velocity/position update portion of the integration.

.. note::

   A Nose-Hoover thermostat will not work well for arbitrary values
   of *Tdamp*\ .  If *Tdamp* is too small, the temperature can fluctuate
   wildly; if it is too large, the temperature will take a very long time
   to equilibrate.  A good choice for many models is a *Tdamp* of around
   100 timesteps.  Note that this is NOT the same as 100 time units for
   most :doc:`units <units>` settings. A simple way to ensure this, is
   via using an :doc:`immediate variable <variable>` expression accessing
   the thermo property 'dt', which is the length of the time step. Example:

.. code-block:: LAMMPS

   fix 1 all nvt temp 300.0 300.0 $(100.0*dt)

----------

The barostat parameters for fix styles *npt* and *nph* is specified
using one or more of the *iso*, *aniso*, *tri*, *x*, *y*, *z*, *xy*,
*xz*, *yz*, and *couple* keywords.  These keywords give you the
ability to specify all 6 components of an external stress tensor, and
to couple various of these components together so that the dimensions
they represent are varied together during a constant-pressure
simulation.

Other barostat-related keywords are *pchain*, *mtk*, *ploop*,
*nreset*, *drag*, and *dilate*, which are discussed below.

Orthogonal simulation boxes have 3 adjustable dimensions (x,y,z).
Triclinic (non-orthogonal) simulation boxes have 6 adjustable
dimensions (x,y,z,xy,xz,yz).  The :doc:`create_box <create_box>`, :doc:`read data <read_data>`, and :doc:`read_restart <read_restart>` commands
specify whether the simulation box is orthogonal or non-orthogonal
(triclinic) and explain the meaning of the xy,xz,yz tilt factors.

The target pressures for each of the 6 components of the stress tensor
can be specified independently via the *x*, *y*, *z*, *xy*, *xz*, *yz*
keywords, which correspond to the 6 simulation box dimensions.  For
each component, the external pressure or tensor component at each
timestep is a ramped value during the run from *Pstart* to *Pstop*\ .
If a target pressure is specified for a component, then the
corresponding box dimension will change during a simulation.  For
example, if the *y* keyword is used, the y-box length will change.  If
the *xy* keyword is used, the xy tilt factor will change.  A box
dimension will not change if that component is not specified, although
you have the option to change that dimension via the :doc:`fix deform <fix_deform>` command.

Note that in order to use the *xy*, *xz*, or *yz* keywords, the
simulation box must be triclinic, even if its initial tilt factors are
0.0.

For all barostat keywords, the *Pdamp* parameter operates like the
*Tdamp* parameter, determining the time scale on which pressure is
relaxed.  For example, a value of 10.0 means to relax the pressure in
a timespan of (roughly) 10 time units (e.g. :math:`\tau` or fs or ps
- see the :doc:`units <units>` command).

.. note::

   A Nose-Hoover barostat will not work well for arbitrary values
   of *Pdamp*\ .  If *Pdamp* is too small, the pressure and volume can
   fluctuate wildly; if it is too large, the pressure will take a very
   long time to equilibrate.  A good choice for many models is a *Pdamp*
   of around 1000 timesteps.  However, note that *Pdamp* is specified in
   time units, and that timesteps are NOT the same as time units for most
   :doc:`units <units>` settings.

The relaxation rate of the barostat is set by its inertia :math:`W`:

.. math::

   W = (N + 1) k_B T_{\rm target} P_{\rm damp}^2

where :math:`N` is the number of atoms, :math:`k_B` is the Boltzmann constant,
and :math:`T_{\rm target}` is the target temperature of the barostat :ref:`(Martyna) <nh-Martyna>`.
If a thermostat is defined, :math:`T_{\rm target}` is the target temperature
of the thermostat. If a thermostat is not defined, :math:`T_{\rm target}`
is set to the current temperature of the system when the barostat is initialized.
If this temperature is too low the simulation will quit with an error.
Note: in previous versions of LAMMPS, :math:`T_{\rm target}` would default to
a value of 1.0 for *lj* units and 300.0 otherwise if the system had a temperature
of exactly zero.

If a thermostat is not specified by this fix, :math:`T_{\rm target}` can be
manually specified using the *Ptemp* parameter. This may be useful if the
barostat is initialized when the current temperature does not reflect the
steady state temperature of the system. This keyword may also be useful in
athermal simulations where the temperature is not well defined.

Regardless of what atoms are in the fix group (the only atoms which
are time integrated), a global pressure or stress tensor is computed
for all atoms.  Similarly, when the size of the simulation box is
changed, all atoms are re-scaled to new positions, unless the keyword
*dilate* is specified with a *dilate-group-ID* for a group that
represents a subset of the atoms.  This can be useful, for example, to
leave the coordinates of atoms in a solid substrate unchanged and
controlling the pressure of a surrounding fluid.  This option should
be used with care, since it can be unphysical to dilate some atoms and
not others, because it can introduce large, instantaneous
displacements between a pair of atoms (one dilated, one not) that are
far from the dilation origin.  Also note that for atoms not in the fix
group, a separate time integration fix like :doc:`fix nve <fix_nve>` or
:doc:`fix nvt <fix_nh>` can be used on them, independent of whether they
are dilated or not.

----------

The *couple* keyword allows two or three of the diagonal components of
the pressure tensor to be "coupled" together.  The value specified
with the keyword determines which are coupled.  For example, *xz*
means the *Pxx* and *Pzz* components of the stress tensor are coupled.
*Xyz* means all 3 diagonal components are coupled.  Coupling means two
things: the instantaneous stress will be computed as an average of the
corresponding diagonal components, and the coupled box dimensions will
be changed together in lockstep, meaning coupled dimensions will be
dilated or contracted by the same percentage every timestep.  The
*Pstart*, *Pstop*, *Pdamp* parameters for any coupled dimensions must
be identical.  *Couple xyz* can be used for a 2d simulation; the *z*
dimension is simply ignored.

----------

The *iso*, *aniso*, and *tri* keywords are simply shortcuts that are
equivalent to specifying several other keywords together.

The keyword *iso* means couple all 3 diagonal components together when
pressure is computed (hydrostatic pressure), and dilate/contract the
dimensions together.  Using "iso Pstart Pstop Pdamp" is the same as
specifying these 4 keywords:

.. parsed-literal::

   x Pstart Pstop Pdamp
   y Pstart Pstop Pdamp
   z Pstart Pstop Pdamp
   couple xyz

The keyword *aniso* means *x*, *y*, and *z* dimensions are controlled
independently using the *Pxx*, *Pyy*, and *Pzz* components of the
stress tensor as the driving forces, and the specified scalar external
pressure.  Using "aniso Pstart Pstop Pdamp" is the same as specifying
these 4 keywords:

.. parsed-literal::

   x Pstart Pstop Pdamp
   y Pstart Pstop Pdamp
   z Pstart Pstop Pdamp
   couple none

The keyword *tri* means *x*, *y*, *z*, *xy*, *xz*, and *yz* dimensions
are controlled independently using their individual stress components
as the driving forces, and the specified scalar pressure as the
external normal stress.  Using "tri Pstart Pstop Pdamp" is the same as
specifying these 7 keywords:

.. parsed-literal::

   x Pstart Pstop Pdamp
   y Pstart Pstop Pdamp
   z Pstart Pstop Pdamp
   xy 0.0 0.0 Pdamp
   yz 0.0 0.0 Pdamp
   xz 0.0 0.0 Pdamp
   couple none

----------

In some cases (e.g. for solids) the pressure (volume) and/or
temperature of the system can oscillate undesirably when a Nose/Hoover
barostat and thermostat is applied.  The optional *drag* keyword will
damp these oscillations, although it alters the Nose/Hoover equations.
A value of 0.0 (no drag) leaves the Nose/Hoover formalism unchanged.
A non-zero value adds a drag term; the larger the value specified, the
greater the damping effect.  Performing a short run and monitoring the
pressure and temperature is the best way to determine if the drag term
is working.  Typically a value between 0.2 to 2.0 is sufficient to
damp oscillations after a few periods. Note that use of the drag
keyword will interfere with energy conservation and will also change
the distribution of positions and velocities so that they do not
correspond to the nominal NVT, NPT, or NPH ensembles.

An alternative way to control initial oscillations is to use chain
thermostats. The keyword *tchain* determines the number of thermostats
in the particle thermostat. A value of 1 corresponds to the original
Nose-Hoover thermostat. The keyword *pchain* specifies the number of
thermostats in the chain thermostatting the barostat degrees of
freedom. A value of 0 corresponds to no thermostatting of the
barostat variables.

The *mtk* keyword controls whether or not the correction terms due to
Martyna, Tuckerman, and Klein are included in the equations of motion
:ref:`(Martyna) <nh-Martyna>`.  Specifying *no* reproduces the original
Hoover barostat, whose volume probability distribution function
differs from the true NPT and NPH ensembles by a factor of 1/V.  Hence
using *yes* is more correct, but in many cases the difference is
negligible.

The keyword *tloop* can be used to improve the accuracy of integration
scheme at little extra cost.  The initial and final updates of the
thermostat variables are broken up into *tloop* sub-steps, each of
length *dt*\ /\ *tloop*\ . This corresponds to using a first-order
Suzuki-Yoshida scheme :ref:`(Tuckerman) <nh-Tuckerman>`.  The keyword *ploop*
does the same thing for the barostat thermostat.

The keyword *nreset* controls how often the reference dimensions used
to define the strain energy are reset.  If this keyword is not used,
or is given a value of zero, then the reference dimensions are set to
those of the initial simulation domain and are never changed. If the
simulation domain changes significantly during the simulation, then
the final average pressure tensor will differ significantly from the
specified values of the external stress tensor.  A value of *nstep*
means that every *nstep* timesteps, the reference dimensions are set
to those of the current simulation domain.

The *scaleyz*, *scalexz*, and *scalexy* keywords control whether or
not the corresponding tilt factors are scaled with the associated box
dimensions when barostatting triclinic periodic cells.  The default
values *yes* will turn on scaling, which corresponds to adjusting the
linear dimensions of the cell while preserving its shape.  Choosing
*no* ensures that the tilt factors are not scaled with the box
dimensions. See below for restrictions and default values in different
situations. In older versions of LAMMPS, scaling of tilt factors was
not performed. The old behavior can be recovered by setting all three
scale keywords to *no*\ .

The *flip* keyword allows the tilt factors for a triclinic box to
exceed half the distance of the parallel box length, as discussed
below.  If the *flip* value is set to *yes*, the bound is enforced by
flipping the box when it is exceeded.  If the *flip* value is set to
*no*, the tilt will continue to change without flipping.  Note that if
applied stress induces large deformations (e.g. in a liquid), this
means the box shape can tilt dramatically and LAMMPS will run less
efficiently, due to the large volume of communication needed to
acquire ghost atoms around a processor's irregular-shaped subdomain.
For extreme values of tilt, LAMMPS may also lose atoms and generate an
error.

The *fixedpoint* keyword specifies the fixed point for barostat volume
changes. By default, it is the center of the box.  Whatever point is
chosen will not move during the simulation.  For example, if the lower
periodic boundaries pass through (0,0,0), and this point is provided
to *fixedpoint*, then the lower periodic boundaries will remain at
(0,0,0), while the upper periodic boundaries will move twice as
far. In all cases, the particle trajectories are unaffected by the
chosen value, except for a time-dependent constant translation of
positions.

If the *update* keyword is used with the *dipole* value, then the
orientation of the dipole moment of each particle is also updated
during the time integration.  This option should be used for models
where a dipole moment is assigned to finite-size particles,
e.g. spheroids via use of the :doc:`atom_style hybrid sphere dipole <atom_style>` command.

The default dipole orientation integrator can be changed to the
Dullweber-Leimkuhler-McLachlan integration scheme
:ref:`(Dullweber) <nh-Dullweber>` when using *update* with the value
*dipole/dlm*\ . This integrator is symplectic and time-reversible,
giving better energy conservation and allows slightly longer timesteps
at only a small additional computational cost.

----------

.. note::

   Using a barostat coupled to tilt dimensions *xy*, *xz*, *yz* can
   sometimes result in arbitrarily large values of the tilt dimensions,
   i.e. a dramatically deformed simulation box.  LAMMPS allows the tilt
   factors to grow a small amount beyond the normal limit of half the box
   length (0.6 times the box length), and then performs a box "flip" to
   an equivalent periodic cell.  See the discussion of the *flip* keyword
   above, to allow this bound to be exceeded, if desired.

The flip operation is described in more detail in the page for
:doc:`fix deform <fix_deform>`.  Both the barostat dynamics and the atom
trajectories are unaffected by this operation.  However, if a tilt
factor is incremented by a large amount (1.5 times the box length) on
a single timestep, LAMMPS can not accommodate this event and will
terminate the simulation with an error. This error typically indicates
that there is something badly wrong with how the simulation was
constructed, such as specifying values of *Pstart* that are too far
from the current stress value, or specifying a timestep that is too
large. Triclinic barostatting should be used with care. This also is
true for other barostat styles, although they tend to be more
forgiving of insults. In particular, it is important to recognize that
equilibrium liquids can not support a shear stress and that
equilibrium solids can not support shear stresses that exceed the
yield stress.

One exception to this rule is if the first dimension in the tilt factor
(x for xy) is non-periodic.  In that case, the limits on the tilt
factor are not enforced, since flipping the box in that dimension does
not change the atom positions due to non-periodicity.  In this mode,
if you tilt the system to extreme angles, the simulation will simply
become inefficient due to the highly skewed simulation box.

.. note::

   Unlike the :doc:`fix temp/berendsen <fix_temp_berendsen>` command
   which performs thermostatting but NO time integration, these fixes
   perform thermostatting/barostatting AND time integration.  Thus you
   should not use any other time integration fix, such as :doc:`fix nve <fix_nve>` on atoms to which this fix is applied.  Likewise,
   fix nvt and fix npt should not normally be used on atoms that also
   have their temperature controlled by another fix - e.g. by :doc:`fix langevin <fix_nh>` or :doc:`fix temp/rescale <fix_temp_rescale>`
   commands.

See the :doc:`Howto thermostat <Howto_thermostat>` and :doc:`Howto barostat <Howto_barostat>` doc pages for a discussion of different
ways to compute temperature and perform thermostatting and
barostatting.

----------

These fixes compute a temperature and pressure each timestep.  To do
this, the thermostat and barostat fixes create their own computes of
style "temp" and "pressure", as if one of these sets of commands had
been issued:

For fix nvt:

.. code-block:: LAMMPS

   compute fix-ID_temp group-ID temp

For fix npt and fix nph:

.. code-block:: LAMMPS

   compute fix-ID_temp all temp
   compute fix-ID_press all pressure fix-ID_temp

For fix nvt, the group for the new temperature compute is the same as
the fix group.  For fix npt and fix nph, the group for both the new
temperature and pressure compute is "all" since pressure is computed
for the entire system.  In the case of fix nph, the temperature
compute is not used for thermostatting, but just for a kinetic-energy
contribution to the pressure.  See the :doc:`compute temp <compute_temp>` and :doc:`compute pressure <compute_pressure>`
commands for details.  Note that the IDs of the new computes are the
fix-ID + underscore + "temp" or fix_ID + underscore + "press".

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo_temp*
and *thermo_press*.  This means you can change the attributes of these
fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command.  Or you can print this
temperature or pressure during thermodynamic output via the
:doc:`thermo_style custom <thermo_style>` command using the appropriate
compute-ID.  It also means that changing attributes of *thermo_temp*
or *thermo_press* will have no effect on this fix.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that remove a "bias" from the
atom velocities.  E.g. to apply the thermostat only to atoms within a
spatial :doc:`region <region>`, or to remove the center-of-mass
velocity from a group of atoms, or to remove the x-component of
velocity from the calculation.

This is not done by default, but only if the :doc:`fix_modify
<fix_modify>` command is used to assign a temperature compute to this
fix that includes such a bias term.  See the doc pages for individual
:doc:`compute temp commands <compute>` to determine which ones include
a bias.  In this case, the thermostat works in the following manner:
bias is removed from each atom, thermostatting is performed on the
remaining thermal degrees of freedom, and the bias is added back in.

----------

These fixes can be used with either the *verlet* or *respa*
:doc:`integrators <run_style>`. When using one of the barostat fixes
with *respa*, LAMMPS uses an integrator constructed
according to the following factorization of the Liouville propagator
(for two rRESPA levels):

.. math::

   \exp \left(\mathrm{i} L \Delta t \right) = & \hat{E}
   \exp \left(\mathrm{i} L_{\rm T\textrm{-}baro} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{\rm T\textrm{-}part} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{\epsilon , 2} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{2}^{(2)} \frac{\Delta t}{2} \right) \\
   &\times \left[
   \exp \left(\mathrm{i} L_{2}^{(1)} \frac{\Delta t}{2n} \right)
   \exp \left(\mathrm{i} L_{\epsilon , 1} \frac{\Delta t}{2n} \right)
   \exp \left(\mathrm{i} L_1 \frac{\Delta t}{n} \right)
   \exp \left(\mathrm{i} L_{\epsilon , 1} \frac{\Delta t}{2n} \right)
   \exp \left(\mathrm{i} L_{2}^{(1)} \frac{\Delta t}{2n} \right)
   \right]^n \\
   &\times
   \exp \left(\mathrm{i} L_{2}^{(2)} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{\epsilon , 2} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{\rm T\textrm{-}part} \frac{\Delta t}{2} \right)
   \exp \left(\mathrm{i} L_{\rm T\textrm{-}baro} \frac{\Delta t}{2} \right) \\
   &+ \mathcal{O} \left(\Delta t^3 \right)

This factorization differs somewhat from that of Tuckerman et al, in
that the barostat is only updated at the outermost rRESPA level,
whereas Tuckerman's factorization requires splitting the pressure into
pieces corresponding to the forces computed at each rRESPA level. In
theory, the latter method will exhibit better numerical stability. In
practice, because Pdamp is normally chosen to be a large multiple of
the outermost rRESPA timestep, the barostat dynamics are not the
limiting factor for numerical stability. Both factorizations are
time-reversible and can be shown to preserve the phase space measure
of the underlying non-Hamiltonian equations of motion.

.. note::

   This implementation has been shown to conserve linear momentum
   up to machine precision under NVT dynamics. Under NPT dynamics,
   for a system with zero initial total linear momentum, the total
   momentum fluctuates close to zero. It may occasionally undergo brief
   excursions to non-negligible values, before returning close to zero.
   Over long simulations, this has the effect of causing the center-of-mass
   to undergo a slow random walk. This can be mitigated by resetting
   the momentum at infrequent intervals using the
   :doc:`fix momentum <fix_momentum>` command.

----------

The fix npt and fix nph commands can be used with rigid bodies or
mixtures of rigid bodies and non-rigid particles (e.g. solvent).  But
there are also :doc:`fix rigid/npt <fix_rigid>` and :doc:`fix rigid/nph <fix_rigid>` commands, which are typically a more natural
choice.  See the page for those commands for more discussion of
the various ways to do this.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes writes the state of all the thermostat and barostat
variables to :doc:`binary restart files <restart>`.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by these fixes.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its thermostatting or barostatting procedure, as described above.
If you do this, note that the kinetic energy derived from the compute
temperature should be consistent with the virial term computed using
all atoms for the pressure.  LAMMPS will warn you if you choose to
compute temperature on a subset of atoms.

.. note::

   If both the *temp* and *press* keywords are used in a single
   thermo_modify command (or in two separate commands), then the order in
   which the keywords are specified is important.  Note that a :doc:`pressure compute <compute_pressure>` defines its own temperature compute as
   an argument when it is specified.  The *temp* keyword will override
   this (for the pressure compute being used by fix npt), but only if the
   *temp* keyword comes after the *press* keyword.  If the *temp* keyword
   comes before the *press* keyword, then the new pressure compute
   specified by the *press* keyword will be unaffected by the *temp*
   setting.

The cumulative energy change in the system imposed by these fixes, via
either thermostatting and/or barostatting, is included in the
:doc:`thermodynamic output <thermo_style>` keywords *ecouple* and
*econserve*.  See the :doc:`thermo_style <thermo_style>` page for
details.

These fixes compute a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the same
cumulative energy change due to this fix described in the previous
paragraph.  The scalar value calculated by this fix is "extensive".

These fixes compute also compute a global vector of quantities, which
can be accessed by various :doc:`output commands <Howto_output>`.  The
vector values are "intensive".

The vector stores internal Nose/Hoover thermostat and barostat
variables.  The number and meaning of the vector values depends on
which fix is used and the settings for keywords *tchain* and *pchain*,
which specify the number of Nose/Hoover chains for the thermostat and
barostat.  If no thermostatting is done, then *tchain* is 0.  If no
barostatting is done, then *pchain* is 0.  In the following list,
"ndof" is 0, 1, 3, or 6, and is the number of degrees of freedom in
the barostat.  Its value is 0 if no barostat is used, else its value
is 6 if any off-diagonal stress tensor component is barostatted, else
its value is 1 if *couple xyz* is used or *couple xy* for a 2d
simulation, otherwise its value is 3.

The order of values in the global vector and their meaning is as
follows.  The notation means there are tchain values for eta, followed
by tchain for eta_dot, followed by ndof for omega, etc:

* eta[tchain] = particle thermostat displacements (unitless)
* eta_dot[tchain] = particle thermostat velocities (1/time units)
* omega[ndof] = barostat displacements (unitless)
* omega_dot[ndof] = barostat velocities (1/time units)
* etap[pchain] = barostat thermostat displacements (unitless)
* etap_dot[pchain] = barostat thermostat velocities (1/time units)
* PE_eta[tchain] = potential energy of each particle thermostat displacement (energy units)
* KE_eta_dot[tchain] = kinetic energy of each particle thermostat velocity (energy units)
* PE_omega[ndof] = potential energy of each barostat displacement (energy units)
* KE_omega_dot[ndof] = kinetic energy of each barostat velocity (energy units)
* PE_etap[pchain] = potential energy of each barostat thermostat displacement (energy units)
* KE_etap_dot[pchain] = kinetic energy of each barostat thermostat velocity (energy units)
* PE_strain[1] = scalar strain energy (energy units)

These fixes can ramp their external temperature and pressure over
multiple runs, using the *start* and *stop* keywords of the
:doc:`run <run>` command.  See the :doc:`run <run>` command for details of
how to do this.

These fixes are not invoked during :doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

*X*, *y*, *z* cannot be barostatted if the associated dimension is not
periodic.  *Xy*, *xz*, and *yz* can only be barostatted if the
simulation domain is triclinic and the second dimension in the keyword
(\ *y* dimension in *xy*\ ) is periodic.  *Z*, *xz*, and *yz*, cannot be
barostatted for 2D simulations.  The :doc:`create_box <create_box>`,
:doc:`read data <read_data>`, and :doc:`read_restart <read_restart>`
commands specify whether the simulation box is orthogonal or
non-orthogonal (triclinic) and explain the meaning of the xy,xz,yz
tilt factors.

For the *temp* keyword, the final Tstop cannot be 0.0 since it would
make the external T = 0.0 at some timestep during the simulation which
is not allowed in the Nose/Hoover formulation.

The *scaleyz yes* and *scalexz yes* keyword/value pairs can not be used
for 2D simulations. *scaleyz yes*, *scalexz yes*, and *scalexy yes* options
can only be used if the second dimension in the keyword is periodic,
and if the tilt factor is not coupled to the barostat via keywords
*tri*, *yz*, *xz*, and *xy*\ .

These fixes can be used with dynamic groups as defined by the
:doc:`group <group>` command.  Likewise they can be used with groups to
which atoms are added or deleted over time, e.g. a deposition
simulation.  However, the conservation properties of the thermostat
and barostat are defined for systems with a static set of atoms.  You
may observe odd behavior if the atoms in a group vary dramatically
over time or the atom count becomes very small.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix_modify <fix_modify>`,
:doc:`run_style <run_style>`

Default
"""""""

The keyword defaults are tchain = 3, pchain = 3, mtk = yes, tloop = 1,
ploop = 1, nreset = 0, drag = 0.0, dilate = all, couple = none,
flip = yes, scaleyz = scalexz = scalexy = yes if periodic in second
dimension and not coupled to barostat, otherwise no.

----------

.. _nh-Martyna:

**(Martyna)** Martyna, Tobias and Klein, J Chem Phys, 101, 4177 (1994).

.. _nh-Parrinello:

**(Parrinello)** Parrinello and Rahman, J Appl Phys, 52, 7182 (1981).

.. _nh-Tuckerman:

**(Tuckerman)** Tuckerman, Alejandre, Lopez-Rendon, Jochim, and
Martyna, J Phys A: Math Gen, 39, 5629 (2006).

.. _nh-Shinoda:

**(Shinoda)** Shinoda, Shiga, and Mikami, Phys Rev B, 69, 134103 (2004).

.. _nh-Dullweber:

**(Dullweber)** Dullweber, Leimkuhler and McLachlan, J Chem Phys, 107,
5840 (1997).
