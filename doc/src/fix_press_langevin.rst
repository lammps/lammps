.. index:: fix press/langevin

fix press/langevin command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID press/langevin keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* press/langevin = style name of this fix command

  .. parsed-literal::

     one or more keyword value pairs may be appended
     keyword = *iso* or *aniso* or *tri* or *x* or *y* or *z* or *xy* or *xz* or *yz* or *couple* or *dilate* or *modulus* or *temp* or *flip*
       *iso* or *aniso* or *tri* values = Pstart Pstop Pdamp
         Pstart,Pstop = scalar external pressure at start/end of run (pressure units)
         Pdamp = pressure damping parameter (time units)
       *x* or *y* or *z* or *xy* or *xz* or *yz* values = Pstart Pstop Pdamp
         Pstart,Pstop = external stress tensor component at start/end of run (pressure units)
         Pdamp = pressure damping parameter
       *flip* value = *yes* or *no* = allow or disallow box flips when it becomes highly skewed
       *couple* = *none* or *xyz* or *xy* or *yz* or *xz*
       *friction* value = Friction coefficient for the barostat (time units)
       *temp* values = Tstart, Tstop, seed
       Tstart, Tstop = target temperature used for the barostat at start/end of run
       seed = seed of the random number generator
       *dilate* value = *all* or *partial*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all press/langevin iso 0.0 0.0 1000.0 temp 300 300 487374
   fix 2 all press/langevin aniso 0.0 0.0 1000.0 temp 100 300 238 dilate partial

Description
"""""""""""

Adjust the pressure of the system by using a Langevin stochastic barostat
:ref:`(Gronbech) <Gronbech>`, which rescales the system volume and
(optionally) the atoms coordinates within the simulation box every
timestep.

The Langevin barostat couple each direction *L* with a pseudo-particle that obeys
the Langevin equation such as:

.. math::

   f_P = & \frac{N k_B T_{target}}{V} + \frac{1}{V d}\sum_{i=1}^{N} \vec r_i \cdot \vec f_i - P_{target} \\
   Q\ddot{L} + \alpha{}\dot{L} = & f_P + \beta(t)\\
   L^{n+1} = & L^{n} + bdt\dot{L}^{n} \frac{bdt^{2}}{2Q} \\
   \dot{L}^{n+1} = & \alpha\dot{L}^{n} + \frac{dt}{2Q}\left(a f^{n}_{P} + f^{n+1}_{P}\right) + \frac{b}{Q}\beta^{n+1} \\
   a = & \frac{1-\frac{\alpha{}dt}{2Q}}{1+\frac{\alpha{}dt}{2Q}} \\
   b = & \frac{1}{1+\frac{\alpha{}dt}{2Q}} \\
   \left< \beta(t)\beta(t') \right> = & 2\alpha k_B Tdt

Where :math:`dt` is the timestep :math:`\dot{L}` and :math:`\ddot{L}` the first
and second derivatives of the coupled direction with regard to time,
:math:`\alpha` is a friction coefficient, :math:`\beta` is a random gaussian
variable and :math:`Q` the effective mass of the coupled pseudoparticle. The
two first terms on the right-hand side of the first equation are the virial
expression of the canonical pressure. It is to be noted that the temperature
used to compute the pressure is not based on the atom velocities but rather on
the canonical
target temperature directly. This temperature is specified using the *temp*
keyword parameter and should be close to the expected target temperature of the
system.

Regardless of what atoms are in the fix group, a global pressure is
computed for all atoms. Similarly, when the size of the simulation
box is changed, all atoms are re-scaled to new positions, unless the
keyword *dilate* is specified with a value of *partial*, in which case
only the atoms in the fix group are re-scaled. The latter can be
useful for leaving the coordinates of atoms in a solid substrate
unchanged and controlling the pressure of a surrounding fluid.

.. note::

   Unlike the :doc:`fix npt <fix_nh>` or :doc:`fix nph <fix_nh>` commands which
   perform Nose-Hoover barostatting AND time integration, this fix does NOT
   perform time integration of the atoms but only of the barostat coupled
   coordinate. It then only modifies the box size and atom coordinates to
   effect barostatting. Thus you must use a separate time integration fix,
   like :doc:`fix nve <fix_nve>` or :doc:`fix nvt <fix_nh>` to actually update
   the positions and velocities of atoms.  This fix can be used in conjunction
   with thermostatting fixes to control the temperature, such as :doc:`fix nvt
   <fix_nh>` or :doc:`fix langevin <fix_langevin>` or :doc:`fix temp/berendsen
   <fix_temp_berendsen>`.

See the :doc:`Howto barostat <Howto_barostat>` page for a
discussion of different ways to perform barostatting.

----------

The barostat is specified using one or more of the *iso*, *aniso*, *tri* *x*,
*y*, *z*, *xy*, *xz*, *yz*, and *couple* keywords.  These keywords give you the
ability to specify the 3 diagonal components of an external stress tensor, and
to couple various of these components together so that the dimensions they
represent are varied together during a constant-pressure simulation.

The target pressures for each of the 6 diagonal components of the stress tensor
can be specified independently via the *x*, *y*, *z*, keywords, which
correspond to the 3 simulation box dimensions, and the *xy*, *xz* and *yz*
keywords which corresponds to the 3 simulation box tilt factors. For each
component, the external pressure or tensor component at each timestep is a
ramped value during the run from *Pstart* to *Pstop*\ . If a target pressure is
specified for a component, then the corresponding box dimension will change
during a simulation.  For example, if the *y* keyword is used, the y-box length
will change.  A box dimension will not change if that component is not
specified, although you have the option to change that dimension via the
:doc:`fix deform <fix_deform>` command.

The *Pdamp* parameter can be seen in the same way as a Nose-Hoover parameter as
it is used to compute the mass of the fictitious particle. Without friction,
the barostat can be compared to a single particle Nose-Hoover barostat and
should follow a similar decay in time. The mass of the barostat is
linked to *Pdamp* by the relation
:math:`Q=(N_{at}+1)\cdot{}k_BT_{target}\cdot{}P_{damp}^2`. Note that *Pdamp*
should be expressed in time units.

.. note::

   As for Berendsen barostat, a Langevin barostat will not work well for
   arbitrary values of *Pdamp*\ .  If *Pdamp* is too small, the pressure and
   volume can fluctuate wildly; if it is too large, the pressure will take a
   very long time to equilibrate.  A good choice for many models is a *Pdamp*
   of around 1000 timesteps.  However, note that *Pdamp* is specified in time
   units, and that timesteps are NOT the same as time units for most
   :doc:`units <units>` settings.

----------

The *temp* keyword sets the temperature to use in the equation of motion of the
barostat. This value is used to compute the value of the force :math:`f_P` in
the equation of motion. It is important to note that this value is not the
instantaneous temperature but a target temperature that ramps from *Tstart* to
*Tstop*. Also the required argument *seed* sets the seed for the random
number generator used in the generation of the random forces.

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

The *iso*, *aniso* and *tri* keywords are simply shortcuts that are
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

The keyword *tri* is the same as *aniso* but also adds the control on the
shear pressure coupled with the tilt factors.

.. parsed-literal::

   x Pstart Pstop Pdamp
   y Pstart Pstop Pdamp
   z Pstart Pstop Pdamp
   xy Pstart Pstop Pdamp
   xz Pstart Pstop Pdamp
   yz Pstart Pstop Pdamp
   couple none

----------

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

----------

The *friction* keyword sets the friction parameter :math:`\alpha` in the
equations of motion of the barostat. For each barostat direction, the value of
:math:`\alpha` depends on both *Pdamp* and *friction*. The value given as a
parameter is the Langevin characteristic time
:math:`\tau_{L}=\frac{Q}{\alpha}` in time units. The langevin time can be understood as a
decorrelation time for the pressure. A long Langevin time value will make the
barostat act as an underdamped oscillator while a short value will make it
act as an overdamped oscillator. The ideal configuration would be to find
the critical parameter of the barostat. Empirically this is observed to
occur for :math:`\tau_{L}\approx{}P_{damp}`.  For this reason, if the *friction*
keyword is not used, the default value *Pdamp* is used for each barostat direction.

----------

This fix computes pressure each timestep. To do
this, the fix creates its own computes of style "pressure",
as if this command had been issued:

.. code-block:: LAMMPS

   compute fix-ID_press group-ID pressure NULL virial

The kinetic contribution to the pressure is taken as the ensemble value
:math:`\frac{Nk_bT}{V}` and computed by the fix itself.

See the :doc:`compute pressure <compute_pressure>` command for details.  Note
that the IDs of the new compute is the fix-ID + underscore + "press" and the
group for the new computes is the same as the fix group.

Note that this is NOT the compute used by thermodynamic output (see the
:doc:`thermo_style <thermo_style>` command) with ID = *thermo_press*. This
means you can change the attributes of this fix's pressure via the
:doc:`compute_modify <compute_modify>` command or print this temperature or
pressure during thermodynamic output via the :doc:`thermo_style custom
<thermo_style>` command using the appropriate compute-ID. It also means that
changing attributes of *thermo_temp* or *thermo_press* will have no effect on
this fix.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *press* option is
supported by this fix.  You can use it to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its pressure calculations.

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.

This fix can ramp its target pressure and temperature over multiple runs, using
the *start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this. It is recommended that
the ramped temperature is the same as the effective temperature of the
thermostatted system. That is, if the system's temperature is ramped by other
commands, it is recommended to do the same with this pressure control.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

Any dimension being adjusted by this fix must be periodic.

Related commands
""""""""""""""""

:doc:`fix press/berendsen <fix_press_berendsen>`,
:doc:`fix nve <fix_nve>`, :doc:`fix nph <fix_nh>`, :doc:`fix npt <fix_nh>`, :doc:`fix langevin <fix_langevin>`,
:doc:`fix_modify <fix_modify>`

Default
"""""""

The keyword defaults are *dilate* = all, *flip* = yes, and *friction* = *Pdamp*.

----------

.. _Gronbech:

**(Gronbech)** Gronbech-Jensen, Farago, J Chem Phys, 141, 194108 (2014).
