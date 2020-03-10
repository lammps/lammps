.. index:: fix nphug

fix nphug command
=================

fix nphug/omp command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nphug keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command

  .. parsed-literal::

     one or more keyword value pairs may be appended
     keyword = *temp* or *iso* or *aniso* or *tri* or *x* or *y* or *z* or *couple* or *tchain* or *pchain* or *mtk* or *tloop* or *ploop* or *nreset* or *drag* or *dilate* or *scaleyz* or *scalexz* or *scalexy*
       *temp* values = Value1 Value2 Tdamp
         Value1, Value2 = Nose-Hoover target temperatures, ignored by Hugoniostat
         Tdamp = temperature damping parameter (time units)
       *iso* or *aniso* or *tri* values = Pstart Pstop Pdamp
         Pstart,Pstop = scalar external pressures, must be equal (pressure units)
         Pdamp = pressure damping parameter (time units)
       *x* or *y* or *z* or *xy* or *yz* or *xz* values = Pstart Pstop Pdamp
         Pstart,Pstop = external stress tensor components, must be equal (pressure units)
         Pdamp = stress damping parameter (time units)
       *couple* = *none* or *xyz* or *xy* or *yz* or *xz*
       *tchain* value = length of thermostat chain (1 = single thermostat)
       *pchain* values = length of thermostat chain on barostat (0 = no thermostat)
       *mtk* value = *yes* or *no* = add in MTK adjustment term or not
       *tloop* value = number of sub-cycles to perform on thermostat
       *ploop* value = number of sub-cycles to perform on barostat thermostat
       *nreset* value = reset reference cell every this many timesteps
       *drag* value = drag factor added to barostat/thermostat (0.0 = no drag)
       *dilate* value = *all* or *partial*
       *scaleyz* value = *yes* or *no* = scale yz with lz
       *scalexz* value = *yes* or *no* = scale xz with lz
       *scalexy* value = *yes* or *no* = scale xy with ly



Examples
""""""""


.. parsed-literal::

   fix myhug all nphug temp 1.0 1.0 10.0 z 40.0 40.0 70.0
   fix myhug all nphug temp 1.0 1.0 10.0 iso 40.0 40.0 70.0 drag 200.0 tchain 1 pchain 0

Description
"""""""""""

This command is a variant of the Nose-Hoover
:doc:`fix npt <fix_nh>` fix style.
It performs time integration of the Hugoniostat equations
of motion developed by Ravelo et al. :ref:`(Ravelo) <Ravelo1>`.
These equations compress the system to a state with average
axial stress or pressure equal to the specified target value
and that satisfies the Rankine-Hugoniot (RH)
jump conditions for steady shocks.

The compression can be performed
either
hydrostatically (using keyword *iso*\ , *aniso*\ , or *tri*\ ) or uniaxially
(using keywords *x*\ , *y*\ , or *z*\ ).  In the hydrostatic case,
the cell dimensions change dynamically so that the average axial stress
in all three directions converges towards the specified target value.
In the uniaxial case, the chosen cell dimension changes dynamically
so that the average
axial stress in that direction converges towards the target value. The
other two cell dimensions are kept fixed (zero lateral strain).

This leads to the following additional restrictions on the keywords:

* One and only one of the following keywords should be used: *iso*\ , *aniso*\ , *tri*\ , *x*\ , *y*\ , *z*
* The specified initial and final target pressures must be the same.
* The keywords *xy*\ , *xz*\ , *yz* may not be used.
* The only admissible value for the couple keyword is *xyz*\ , which has the same effect as keyword *iso*
* The *temp* keyword must be used to specify the time constant for kinetic energy relaxation, but initial and final target temperature values are ignored.

Essentially, a Hugoniostat simulation is an NPT simulation in which the
user-specified target temperature is replaced with a time-dependent
target temperature Tt obtained from the following equation:

.. math::

   T_t - T = \frac{\left(\frac{1}{2}\left(P + P_0\right)\left(V_0 - V\right) + E_0 - E\right)}{N_{dof} k_B } = \Delta


where *T* and :math:`T_t` are the instantaneous and target temperatures,
*P* and :math:`P_0` are the instantaneous and reference pressures or axial stresses,
depending on whether hydrostatic or uniaxial compression is being
performed, *V* and :math:`V_0` are the instantaneous and reference volumes,
*E* and :math:`E_0` are the instantaneous and reference internal energy (potential
plus kinetic), :math:`N_{dof}` is the number of degrees of freedom used in the
definition of temperature, and :math:`k_B` is the Boltzmann constant. :math:`\Delta` is the
negative deviation of the instantaneous temperature from the target temperature.
When the system reaches a stable equilibrium, the value of :math:`\Delta` should
fluctuate about zero.

The values of :math:`E_0`, :math:`V_0`, and :math:`P_0` are the instantaneous values at the start of
the simulation. These can be overridden using the fix\_modify keywords *e0*\ ,
*v0*\ , and *p0* described below.


----------


.. note::

   Unlike the :doc:`fix temp/berendsen <fix_temp_berendsen>` command
   which performs thermostatting but NO time integration, this fix
   performs thermostatting/barostatting AND time integration.  Thus you
   should not use any other time integration fix, such as :doc:`fix nve <fix_nve>` on atoms to which this fix is applied.  Likewise,
   this fix should not be used on atoms that have their temperature
   controlled by another fix - e.g. by :doc:`fix langevin <fix_nh>` or :doc:`fix temp/rescale <fix_temp_rescale>` commands.


----------


This fix computes a temperature and pressure at each timestep.  To do
this, the fix creates its own computes of style "temp" and "pressure",
as if one of these two sets of commands had been issued:


.. parsed-literal::

   compute fix-ID_temp group-ID temp
   compute fix-ID_press group-ID pressure fix-ID_temp

   compute fix-ID_temp all temp
   compute fix-ID_press all pressure fix-ID_temp

See the :doc:`compute temp <compute_temp>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + underscore + "temp" or fix\_ID
+ underscore + "press".  The group for
the new computes is "all" since pressure is computed for the entire
system.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo\_temp*
and *thermo\_press*.  This means you can change the attributes of this
fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command or print this temperature
or pressure during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* or
*thermo\_press* will have no effect on this fix.


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


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the values of :math:`E_0`, :math:`V_0`, and :math:`P_0`,
as well as the state of all the thermostat and barostat variables to
:doc:`binary restart files <restart>`.  See the :doc:`read_restart
<read_restart>` command for info on how to re-specify a fix in an input
script that reads a restart file, so that the operation of the fix
continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *e0*\ , *v0* and *p0* keywords can be
used to define the values of :math:`E_0`, :math:`V_0`, and
:math:`P_0`. Note the the values for *e0* and *v0* are extensive, and so
must correspond to the total energy and volume of the entire system, not
energy and volume per atom. If any of these quantities are not
specified, then the instantaneous value in the system at the start of
the simulation is used.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by these fixes.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its thermostatting or barostatting procedure, as described above.
If you do this, note that the kinetic energy derived from the compute
temperature should be consistent with the virial term computed using
all atoms for the pressure.  LAMMPS will warn you if you choose to
compute temperature on a subset of atoms.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by these
fixes to add the energy change induced by Nose/Hoover thermostatting
and barostatting to the system's potential energy as part of
:doc:`thermodynamic output <thermo_style>`. Either way, this energy is \*not\*
included in the definition of internal energy E when calculating the value
of Delta in the above equation.

These fixes compute a global scalar and a global vector of quantities,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar value calculated by these fixes is "extensive"; the vector
values are "intensive".

The scalar is the cumulative energy change due to the fix.

The vector stores three quantities unique to this fix (:math:`\Delta`, Us, and up),
followed by all the internal Nose/Hoover thermostat and barostat
variables defined for :doc:`fix npt <fix_nh>`. Delta is the deviation
of the temperature from the target temperature, given by the above equation.
Us and up are the shock and particle velocity corresponding to a steady
shock calculated from the RH conditions. They have units of distance/time.

Restrictions
""""""""""""


This fix style is part of the SHOCK package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

All the usual restrictions for :doc:`fix npt <fix_nh>` apply,
plus the additional ones mentioned above.

Related commands
""""""""""""""""

:doc:`fix msst <fix_msst>`, :doc:`fix npt <fix_nh>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

The keyword defaults are the same as those for :doc:`fix npt <fix_nh>`


----------


.. _Ravelo1:



**(Ravelo)** Ravelo, Holian, Germann and Lomdahl, Phys Rev B, 70, 014103 (2004).
