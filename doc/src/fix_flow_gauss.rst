.. index:: fix flow/gauss

fix flow/gauss command
======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID flow/gauss xflag yflag zflag keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* flow/gauss = style name of this fix command
* xflag,yflag,zflag = 0 or 1
  
  .. parsed-literal::
  
         0 = do not conserve current in this dimension
         1 = conserve current in this dimension

* zero or more keyword/value pairs may be appended
* keyword = *energy*
  
  .. parsed-literal::
  
       *energy* value = no or yes
         no = do not compute work done by this fix
         yes = compute work done by this fix



Examples
""""""""


.. parsed-literal::

   fix GD fluid flow/gauss 1 0 0
   fix GD fluid flow/gauss 1 1 1 energy yes

Description
"""""""""""

This fix implements the Gaussian dynamics (GD) method to simulate a
system at constant mass flux :ref:`(Strong) <Strong>`. GD is a
nonequilibrium molecular dynamics simulation method that can be used
to study fluid flows through pores, pipes, and channels. In its
original implementation GD was used to compute the pressure required
to achieve a fixed mass flux through an opening.  The flux can be
conserved in any combination of the directions, x, y, or z, using
xflag,yflag,zflag. This fix does not initialize a net flux through a
system, it only conserves the center-of-mass momentum that is present
when the fix is declared in the input script. Use the
:doc:`velocity <velocity>` command to generate an initial center-of-mass
momentum.

GD applies an external fluctuating gravitational field that acts as a
driving force to keep the system away from equilibrium. To maintain
steady state, a profile-unbiased thermostat must be implemented to
dissipate the heat that is added by the driving force. :doc:`Compute temp/profile <compute_temp_profile>` can be used to implement a
profile-unbiased thermostat.

A common use of this fix is to compute a pressure drop across a pipe,
pore, or membrane. The pressure profile can be computed in LAMMPS with
:doc:`compute stress/atom <compute_stress_atom>` and :doc:`fix ave/chunk <fix_ave_chunk>`, or with the hardy method in :doc:`fix atc <fix_atc>`. Note that the simple :doc:`compute stress/atom <compute_stress_atom>` method is only accurate away
from inhomogeneities in the fluid, such as fixed wall atoms. Further,
the computed pressure profile must be corrected for the acceleration
applied by GD before computing a pressure drop or comparing it to
other methods, such as the pump method :ref:`(Zhu) <Zhu>`. The pressure
correction is discussed and described in :ref:`(Strong) <Strong>`.

For a complete example including the considerations discussed
above, see the examples/USER/flow\_gauss directory.

.. note::

   Only the flux of the atoms in group-ID will be conserved. If the
   velocities of the group-ID atoms are coupled to the velocities of
   other atoms in the simulation, the flux will not be conserved. For
   example, in a simulation with fluid atoms and harmonically constrained
   wall atoms, if a single thermostat is applied to group *all*\ , the
   fluid atom velocities will be coupled to the wall atom velocities, and
   the flux will not be conserved. This issue can be avoided by
   thermostatting the fluid and wall groups separately.

Adding an acceleration to atoms does work on the system. This added
energy can be optionally subtracted from the potential energy for the
thermodynamic output (see below) to check that the timestep is small
enough to conserve energy. Since the applied acceleration is
fluctuating in time, the work cannot be computed from a potential. As
a result, computing the work is slightly more computationally
expensive than usual, so it is not performed by default. To invoke the
work calculation, use the *energy* keyword. The
:doc:`fix_modify <fix_modify>` *energy* option also invokes the work
calculation, and overrides an *energy no* setting here. If neither
*energy yes* or *fix\_modify energy yes* are set, the global scalar
computed by the fix will return zero.

.. note::

   In order to check energy conservation, any other fixes that do
   work on the system must have *fix\_modify energy yes* set as well. This
   includes thermostat fixes and any constraints that hold the positions
   of wall atoms fixed, such as :doc:`fix spring/self <fix_spring_self>`.

If this fix is used in a simulation with the :doc:`rRESPA <run_style>`
integrator, the applied acceleration must be computed and applied at the same
rRESPA level as the interactions between the flowing fluid and the obstacle.
The rRESPA level at which the acceleration is applied can be changed using
the :doc:`fix_modify <fix_modify>` *respa* option discussed below. If the
flowing fluid and the obstacle interact through multiple interactions that are
computed at different rRESPA levels, then there must be a separate flow/gauss
fix for each level. For example, if the flowing fluid and obstacle interact
through pairwise and long-range Coulomb interactions, which are computed at
rRESPA levels 3 and 4, respectively, then there must be two separate
flow/gauss fixes, one that specifies *fix\_modify respa 3* and one with
*fix\_modify respa 4*.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to subtract the work done from the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows the user to set at which level of the :doc:`rRESPA <run_style>`
integrator the fix computes and adds the external acceleration. Default is the
outermost level.

This fix computes a global scalar and a global 3-vector of forces,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar is the negative of the work done on the system, see above
discussion.  The vector is the total force that this fix applied to
the group of atoms on the current timestep.  The scalar and vector
values calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix addforce <fix_addforce>`, :doc:`compute temp/profile <compute_temp_profile>`, :doc:`velocity <velocity>`

Default
"""""""

The option default for the *energy* keyword is energy = no.


----------


.. _Strong:



**(Strong)** Strong and Eaves, J. Phys. Chem. B 121, 189 (2017).

.. _Evans2:



**(Evans)** Evans and Morriss, Phys. Rev. Lett. 56, 2172 (1986).

.. _Zhu:



**(Zhu)** Zhu, Tajkhorshid, and Schulten, Biophys. J. 83, 154 (2002).


