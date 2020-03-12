.. index:: fix press/berendsen

fix press/berendsen command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID press/berendsen keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* press/berendsen = style name of this fix command

  .. parsed-literal::

     one or more keyword value pairs may be appended
     keyword = *iso* or *aniso* or *x* or *y* or *z* or *couple* or *dilate* or *modulus*
       *iso* or *aniso* values = Pstart Pstop Pdamp
         Pstart,Pstop = scalar external pressure at start/end of run (pressure units)
         Pdamp = pressure damping parameter (time units)
       *x* or *y* or *z* values = Pstart Pstop Pdamp
         Pstart,Pstop = external stress tensor component at start/end of run (pressure units)
         Pdamp = stress damping parameter (time units)
       *couple* = *none* or *xyz* or *xy* or *yz* or *xz*
       *modulus* value = bulk modulus of system (pressure units)
       *dilate* value = *all* or *partial*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all press/berendsen iso 0.0 0.0 1000.0
   fix 2 all press/berendsen aniso 0.0 0.0 1000.0 dilate partial

Description
"""""""""""

Reset the pressure of the system by using a Berendsen barostat
:ref:`(Berendsen) <Berendsen1>`, which rescales the system volume and
(optionally) the atoms coordinates within the simulation box every
timestep.

Regardless of what atoms are in the fix group, a global pressure is
computed for all atoms.  Similarly, when the size of the simulation
box is changed, all atoms are re-scaled to new positions, unless the
keyword *dilate* is specified with a value of *partial*\ , in which case
only the atoms in the fix group are re-scaled.  The latter can be
useful for leaving the coordinates of atoms in a solid substrate
unchanged and controlling the pressure of a surrounding fluid.

.. note::

   Unlike the :doc:`fix npt <fix_nh>` or :doc:`fix nph <fix_nh>`
   commands which perform Nose/Hoover barostatting AND time integration,
   this fix does NOT perform time integration.  It only modifies the box
   size and atom coordinates to effect barostatting.  Thus you must use a
   separate time integration fix, like :doc:`fix nve <fix_nve>` or :doc:`fix nvt <fix_nh>` to actually update the positions and velocities of
   atoms.  This fix can be used in conjunction with thermostatting fixes
   to control the temperature, such as :doc:`fix nvt <fix_nh>` or :doc:`fix langevin <fix_langevin>` or :doc:`fix temp/berendsen <fix_temp_berendsen>`.

See the :doc:`Howto baroostat <Howto_barostat>` doc page for a
discussion of different ways to perform barostatting.

----------

The barostat is specified using one or more of the *iso*\ , *aniso*\ ,
*x*\ , *y*\ , *z*\ , and *couple* keywords.  These keywords give you the
ability to specify the 3 diagonal components of an external stress
tensor, and to couple various of these components together so that the
dimensions they represent are varied together during a
constant-pressure simulation.  Unlike the :doc:`fix npt <fix_nh>` and
:doc:`fix nph <fix_nh>` commands, this fix cannot be used with triclinic
(non-orthogonal) simulation boxes to control all 6 components of the
general pressure tensor.

The target pressures for each of the 3 diagonal components of the
stress tensor can be specified independently via the *x*\ , *y*\ , *z*\ ,
keywords, which correspond to the 3 simulation box dimensions.  For
each component, the external pressure or tensor component at each
timestep is a ramped value during the run from *Pstart* to *Pstop*\ .
If a target pressure is specified for a component, then the
corresponding box dimension will change during a simulation.  For
example, if the *y* keyword is used, the y-box length will change.  A
box dimension will not change if that component is not specified,
although you have the option to change that dimension via the :doc:`fix deform <fix_deform>` command.

For all barostat keywords, the *Pdamp* parameter determines the time
scale on which pressure is relaxed.  For example, a value of 10.0
means to relax the pressure in a timespan of (roughly) 10 time units
(tau or fmsec or psec - see the :doc:`units <units>` command).

.. note::

   A Berendsen barostat will not work well for arbitrary values of
   *Pdamp*\ .  If *Pdamp* is too small, the pressure and volume can
   fluctuate wildly; if it is too large, the pressure will take a very
   long time to equilibrate.  A good choice for many models is a *Pdamp*
   of around 1000 timesteps.  However, note that *Pdamp* is specified in
   time units, and that timesteps are NOT the same as time units for most
   :doc:`units <units>` settings.

.. note::

   The relaxation time is actually also a function of the bulk
   modulus of the system (inverse of isothermal compressibility).  The
   bulk modulus has units of pressure and is the amount of pressure that
   would need to be applied (isotropically) to reduce the volume of the
   system by a factor of 2 (assuming the bulk modulus was a constant,
   independent of density, which it's not).  The bulk modulus can be set
   via the keyword *modulus*\ .  The *Pdamp* parameter is effectively
   multiplied by the bulk modulus, so if the pressure is relaxing faster
   than expected or desired, increasing the bulk modulus has the same
   effect as increasing *Pdamp*\ .  The converse is also true.  LAMMPS does
   not attempt to guess a correct value of the bulk modulus; it just uses
   10.0 as a default value which gives reasonable relaxation for a
   Lennard-Jones liquid, but will be way off for other materials and way
   too small for solids.  Thus you should experiment to find appropriate
   values of *Pdamp* and/or the *modulus* when using this fix.

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
*Pstart*\ , *Pstop*\ , *Pdamp* parameters for any coupled dimensions must
be identical.  *Couple xyz* can be used for a 2d simulation; the *z*
dimension is simply ignored.

----------

The *iso* and *aniso* keywords are simply shortcuts that are
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

The keyword *aniso* means *x*\ , *y*\ , and *z* dimensions are controlled
independently using the *Pxx*\ , *Pyy*\ , and *Pzz* components of the
stress tensor as the driving forces, and the specified scalar external
pressure.  Using "aniso Pstart Pstop Pdamp" is the same as specifying
these 4 keywords:

.. parsed-literal::

   x Pstart Pstop Pdamp
   y Pstart Pstop Pdamp
   z Pstart Pstop Pdamp
   couple none

----------

This fix computes a temperature and pressure each timestep.  To do
this, the fix creates its own computes of style "temp" and "pressure",
as if these commands had been issued:

.. code-block:: LAMMPS

   compute fix-ID_temp group-ID temp
   compute fix-ID_press group-ID pressure fix-ID_temp

See the :doc:`compute temp <compute_temp>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + underscore + "temp" or fix\_ID
+ underscore + "press", and the group for the new computes is the same
as the fix group.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo\_temp*
and *thermo\_press*.  This means you can change the attributes of this
fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command or print this temperature
or pressure during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* or
*thermo\_press* will have no effect on this fix.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by this fix.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its temperature and pressure calculations.  If you do this, note
that the kinetic energy derived from the compute temperature should be
consistent with the virial term computed using all atoms for the
pressure.  LAMMPS will warn you if you choose to compute temperature
on a subset of atoms.

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.

This fix can ramp its target pressure over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

Any dimension being adjusted by this fix must be periodic.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nph <fix_nh>`, :doc:`fix npt <fix_nh>`, :doc:`fix temp/berendsen <fix_temp_berendsen>`,
:doc:`fix_modify <fix_modify>`

Default
"""""""

The keyword defaults are dilate = all, modulus = 10.0 in units of
pressure for whatever :doc:`units <units>` are defined.

----------

.. _Berendsen1:

**(Berendsen)** Berendsen, Postma, van Gunsteren, DiNola, Haak, J Chem
Phys, 81, 3684 (1984).
