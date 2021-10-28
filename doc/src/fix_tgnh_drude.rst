.. index:: fix tgnvt/drude
.. index:: fix tgnpt/drude

fix tgnvt/drude command
=======================

fix tgnpt/drude command
=======================


Syntax
""""""

.. parsed-literal::

   fix ID group-ID style_name keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *tgnvt/drude* or *tgnpt/drude*
* one or more keyword/values pairs may be appended

  .. parsed-literal::

     keyword = *temp* *iso* or *aniso* or *tri* or *x* or *y* or *z* or *xy* or *yz* or *xz* or *couple* or *tchain* or *pchain* or *mtk* or *tloop* or *ploop* or *nreset* or *scalexy* or *scaleyz* or *scalexz* or *flip* or *fixedpoint*
       *temp* values = Tstart Tstop Tdamp Tdrude Tdamp_drude
         Tstart, Tstop = external temperature at start/end of run (temperature units)
         Tdamp = temperature damping parameter (time units)
         Tdrude = desired temperature of Drude oscillators (temperature units)
         Tdamp_drude = temperature damping parameter for Drude oscillators (time units)
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
       *scalexy* value = *yes* or *no* = scale xy with ly
       *scaleyz* value = *yes* or *no* = scale yz with lz
       *scalexz* value = *yes* or *no* = scale xz with lz
       *flip* value = *yes* or *no* = allow or disallow box flips when it becomes highly skewed
       *fixedpoint* values = x y z
         x,y,z = perform barostat dilation/contraction around this point (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   comm_modify vel yes
   fix 1 all tgnvt/drude temp 300.0 300.0 100.0 1.0 20.0
   fix 1 water tgnpt/drude temp 300.0 300.0 100.0 1.0 20.0 iso 0.0 0.0 1000.0
   fix 2 jello tgnpt/drude temp 300.0 300.0 100.0 1.0 20.0 tri 5.0 5.0 1000.0
   fix 2 ice tgnpt/drude temp 250.0 250.0 100.0 1.0 20.0 x 1.0 1.0 0.5 y 2.0 2.0 0.5 z 3.0 3.0 0.5 yz 0.1 0.1 0.5 xz 0.2 0.2 0.5 xy 0.3 0.3 0.5 nreset 1000

Example input scripts available: examples/PACKAGES/drude

Description
"""""""""""

These commands are variants of the Nose-Hoover fix styles :doc:`fix nvt
<fix_nh>` and :doc:`fix npt <fix_nh>` for thermalized Drude polarizable
models.  They apply temperature-grouped Nose-Hoover thermostat (TGNH)
proposed by :ref:`(Son) <tgnh-Son>`.  When there are fast vibrational
modes with frequencies close to Drude oscillators (e.g. double bonds or
out-of-plane torsions), this thermostat can provide better kinetic
energy equipartitioning.

The difference between TGNH and the original Nose-Hoover thermostat is that,
TGNH separates the kinetic energy of the group into three contributions:
molecular center of mass (COM) motion,
motion of COM of atom-Drude pairs or non-polarizable atoms relative to molecular COM,
and relative motion of atom-Drude pairs.
An independent Nose-Hoover chain is applied to each type of motion.
The temperatures for these three types of motion are denoted as
molecular translational temperature (:math:`T_\mathrm{M}`), real atomic temperature (:math:`T_\mathrm{R}`) and Drude temperature (:math:`T_\mathrm{D}`),
which are defined in terms of their associated degrees of freedom (DOF):

.. math::

    T_\mathrm{M}=\frac{\Sigma_{i}^{N_\mathrm{mol}} M_i V_i^2}{3 \left ( N_\mathrm{mol} - \frac{N_\mathrm{mol}}{N_\mathrm{mol,sys}} \right ) k_\mathrm{B}}

.. math::

    T_\mathrm{R}=\frac{\Sigma_{i}^{N_\mathrm{real}} m_i (v_i-v_{M,i})^2}{(N_\mathrm{DOF} - 3 N_\mathrm{mol} + 3 \frac{N_\mathrm{mol}}{N_\mathrm{mol,sys}} - 3 N_\mathrm{drude}) k_\mathrm{B}}

.. math::

    T_\mathrm{D}=\frac{\Sigma_{i}^{N_\mathrm{drude}} m_i^{\prime} v_i^{\prime 2}}{3 N_\mathrm{drude} k_\mathrm{B}}

Here :math:`N_\mathrm{mol}` and :math:`N_\mathrm{mol,sys}` are the numbers of molecules in the group and in the whole system, respectively.
:math:`N_\mathrm{real}` is the number of atom-Drude pairs and non-polarizable atoms in the group.
:math:`N_\mathrm{drude}` is the number of Drude particles in the group.
:math:`N_\mathrm{DOF}` is the DOF of the group.
:math:`M_i` and :math:`V_i` are the mass and the COM velocity of the i-th molecule.
:math:`m_i` is the mass of the i-th atom-Drude pair or non-polarizable atom.
:math:`v_i` is the velocity of COM of i-th atom-Drude pair or non-polarizable atom.
:math:`v_{M,i}` is the COM velocity of the molecule the i-th atom-Drude pair or non-polarizable atom belongs to.
:math:`m_i^\prime` and :math:`v_i^\prime` are the reduced mass and the relative velocity of the i-th atom-Drude pair.

.. note::

   These fixes require that each atom knows whether it is a Drude particle or
   not.  You must therefore use the :doc:`fix drude <fix_drude>` command to
   specify the Drude status of each atom type.

   Because the TGNH thermostat thermostats the molecular COM motion,
   all atoms belonging to the same molecule must be in the same group.
   That is, these fixes can not be applied to a subset of a molecule.

   For this fix to act correctly, ghost atoms need to know their velocity.
   You must use the :doc:`comm_modify <comm_modify>` command to enable this.

   These fixes assume that the translational DOF of the whole system is removed.
   It is therefore recommended to invoke :doc:`fix momentum <fix_momentum>` command so that the :math:`T_\mathrm{M}` is calculated correctly.

----------

The thermostat parameters are specified using the *temp* keyword.
The thermostat is applied to only the translational DOF
for the particles.  The translational DOF can also have
a bias velocity removed before thermostatting takes place; see the
description below.  The desired temperature for molecular and real atomic motion is a
ramped value during the run from *Tstart* to *Tstop*\ .  The *Tdamp*
parameter is specified in time units and determines how rapidly the
temperature is relaxed.  For example, a value of 10.0 means to relax
the temperature in a timespan of (roughly) 10 time units (e.g. :math:`\tau`
or fs or ps - see the :doc:`units <units>` command).
The parameter *Tdrude* is the desired temperature for Drude motion at each timestep.
Similar to *Tdamp*, the *Tdamp_drude* parameter determines the relaxation speed for Drude motion.
Fix group are the only ones whose velocities and positions are updated
by the velocity/position update portion of the integration.
Other thermostat-related keywords are *tchain*\  and *tloop*,
which are detailed in :doc:`fix nvt <fix_nh>`.

.. note::

   A Nose-Hoover thermostat will not work well for arbitrary values
   of *Tdamp*\ .  If *Tdamp* is too small, the temperature can fluctuate
   wildly; if it is too large, the temperature will take a very long time
   to equilibrate.  A good choice for many models is a *Tdamp* of around
   100 timesteps.  A smaller *Tdamp_drude* value would be required
   to maintain Drude motion at low temperature.

.. code-block:: LAMMPS

   fix 1 all nvt temp 300.0 300.0 $(100.0*dt) 1.0 $(20.0*dt)

----------

The barostat parameters for fix style *tgnpt/drude* is specified
using one or more of the *iso*, *aniso*, *tri*, *x*, *y*, *z*, *xy*,
*xz*, *yz*, and *couple* keywords.  These keywords give you the
ability to specify all 6 components of an external stress tensor, and
to couple various of these components together so that the dimensions
they represent are varied together during a constant-pressure
simulation. Other barostat-related keywords are *pchain*, *mtk*, *ploop*,
*nreset*, *scalexy*, *scaleyz*, *scalexz*, *flip*\ and *fixedpoint*.
The meaning of barostat parameters are detailed in :doc:`fix npt <fix_nh>`.

Regardless of what atoms are in the fix group (the only atoms which
are time integrated), a global pressure or stress tensor is computed
for all atoms.  Similarly, when the size of the simulation box is
changed, all atoms are re-scaled to new positions.

.. note::

   Unlike the :doc:`fix temp/berendsen <fix_temp_berendsen>` command
   which performs thermostatting but NO time integration, these fixes
   perform thermostatting/barostatting AND time integration.  Thus you
   should not use any other time integration fix, such as :doc:`fix nve <fix_nve>` on atoms to which this fix is applied.
   Likewise, these fixes should not be used on atoms that also
   have their temperature controlled by another fix - e.g. by :doc:`fix langevin/drude <fix_langevin_drude>` command.

See the :doc:`Howto thermostat <Howto_thermostat>` and :doc:`Howto barostat <Howto_barostat>` doc pages for a discussion of different
ways to compute temperature and perform thermostatting and
barostatting.

----------

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

.. note::

   However, not all temperature compute commands are valid to be used
   with these fixes.  Precisely, only temperature compute that does
   not modify the DOF of the group can be used.  E.g. :doc:`compute
   temp/ramp <compute_temp_ramp>` and :doc:`compute viscosity/cos
   <compute_viscosity_cos>` compute the kinetic energy after remove a
   velocity gradient without affecting the DOF of the group, then they
   can be invoked in this way.  In contrast, :doc:`compute
   temp/partial <compute_temp_partial>` may remove the DOF at one or
   more dimensions, therefore it cannot be used with these fixes.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes writes the state of all the thermostat and barostat
variables to :doc:`binary restart files <restart>`.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by these fixes.  You can use them to assign a :doc:`compute
<compute>` you have defined to this fix which will be used in its
thermostatting or barostatting procedure, as described above.  If you
do this, note that the kinetic energy derived from the compute
temperature should be consistent with the virial term computed using
all atoms for the pressure.  LAMMPS will warn you if you choose to
compute temperature on a subset of atoms.

.. note::

   If both the *temp* and *press* keywords are used in a single
   thermo_modify command (or in two separate commands), then the order
   in which the keywords are specified is important.  Note that a
   :doc:`pressure compute <compute_pressure>` defines its own
   temperature compute as an argument when it is specified.  The
   *temp* keyword will override this (for the pressure compute being
   used by fix npt), but only if the *temp* keyword comes after the
   *press* keyword.  If the *temp* keyword comes before the *press*
   keyword, then the new pressure compute specified by the *press*
   keyword will be unaffected by the *temp* setting.

The cumulative energy change in the system imposed by these fixes, due
to thermostatting and/or barostatting, are included in the
:doc:`thermodynamic output <thermo_style>` keywords *ecouple* and
*econserve*.  See the :doc:`thermo_style <thermo_style>` page for
details.

These fixes compute a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the same
cumulative energy change due to this fix described in the previous
paragraph.  The scalar value calculated by this fix is "extensive".

These fixes also compute a global vector of quantities, which can be
accessed by various :doc:`output commands <Howto_output>`.  The vector
values are "intensive".  The vector stores the three temperatures
:math:`T_\mathrm{M}`, :math:`T_\mathrm{R}` and :math:`T_\mathrm{D}`.

These fixes can ramp their external temperature and pressure over
multiple runs, using the *start* and *stop* keywords of the :doc:`run
<run>` command.  See the :doc:`run <run>` command for details of how
to do this.

These fixes are not invoked during :doc:`energy minimization
<minimize>`.

----------

Restrictions
""""""""""""

These fixes are only available when LAMMPS was built with the
DRUDE package.  These fixes cannot be used with dynamic groups as
defined by the :doc:`group <group>` command.  These fixes cannot be
used in 2D simulations.

*X*, *y*, *z* cannot be barostatted if the associated dimension is not
periodic.  *Xy*, *xz*, and *yz* can only be barostatted if the
simulation domain is triclinic and the second dimension in the keyword
(\ *y* dimension in *xy*\ ) is periodic.  The :doc:`create_box <create_box>`,
:doc:`read data <read_data>`, and :doc:`read_restart <read_restart>`
commands specify whether the simulation box is orthogonal or
non-orthogonal (triclinic) and explain the meaning of the xy,xz,yz
tilt factors.

For the *temp* keyword, the final *Tstop* cannot be 0.0 since it would
make the external T = 0.0 at some timestep during the simulation which
is not allowed in the Nose/Hoover formulation.

The *scaleyz yes*, *scalexz yes*, and *scalexy yes* options
can only be used if the second dimension in the keyword is periodic,
and if the tilt factor is not coupled to the barostat via keywords
*tri*, *yz*, *xz*, and *xy*\ .

Related commands
""""""""""""""""

:doc:`fix drude <fix_drude>`, :doc:`fix nvt <fix_nh>`, :doc:`fix_npt <fix_nh>`,
:doc:`fix_modify <fix_modify>`

Default
"""""""

The keyword defaults are tchain = 3, pchain = 3, mtk = yes, tloop = 1,
ploop = 1, nreset = 0, couple = none,
flip = yes, scaleyz = scalexz = scalexy = yes if periodic in second
dimension and not coupled to barostat, otherwise no.

----------

.. _tgnh-Son:

**(Son)** Son, McDaniel, Cui and Yethiraj, J Phys Chem Lett, 10, 7523 (2019).
