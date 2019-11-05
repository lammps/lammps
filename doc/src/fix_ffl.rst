.. index:: fix ffl

fix ffl command
===============

Syntax
""""""


.. parsed-literal::

   fix ID id-group ffl tau Tstart Tstop seed [flip-type]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ffl = style name of this fix command
* tau = thermostat parameter (positive real)
* Tstart, Tstop = temperature ramp during the run
* seed = random number seed to use for generating noise (positive integer)
* one more value may be appended
  
  .. parsed-literal::
  
         flip-type  = determines the flipping type, can be chosen between rescale - no_flip - hard - soft, if no flip type is given, rescale will be chosen by default



Examples
""""""""


.. parsed-literal::

   fix 3 boundary ffl 10 300 300 31415
   fix 1 all ffl 100 500 500 9265 soft

Description
"""""""""""

Apply a Fast-Forward Langevin Equation (FFL) thermostat as described
in :ref:`(Hijazi) <Hijazi>`. Contrary to
:doc:`fix langevin <fix_langevin>`, this fix performs both
thermostatting and evolution of the Hamiltonian equations of motion, so it
should not be used together with :doc:`fix nve <fix_nve>` -- at least not
on the same atom groups.

The time-evolution of a single particle undergoing Langevin dynamics is described
by the equations


.. math::

   \begin{equation} \frac {dq}{dt} = \frac{p}{m}, \end{equation}


.. math::

   \begin{equation} \frac {dp}{dt} = -\gamma p + W + F, \end{equation}

where :math:`F` is the physical force, :math:`\gamma` is the friction coefficient, and :math:`W` is a
Gaussian random force.

The friction coefficient is the inverse of the thermostat parameter : :math:`\gamma = 1/\tau`, with :math:`\tau` the thermostat parameter *tau*\ .
The thermostat parameter is given in the time units, :math:`\gamma` is in inverse time units.

Equilibrium sampling a temperature T is obtained by specifying the
target value as the *Tstart* and *Tstop* arguments, so that the internal
constants depending on the temperature are computed automatically.

The random number *seed* must be a positive integer.  A Marsaglia random
number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.

The flipping type *flip-type* can be chosen between 4 types described in
:ref:`(Hijazi) <Hijazi>`. The flipping operation occurs during the thermostatting
step and it flips the momenta of the atoms. If no\_flip is chosen, no flip
will be executed and the integration will be the same as a standard
Langevin thermostat :ref:`(Bussi) <Bussi3>`. The other flipping types are : rescale - hard - soft.

**Restart, fix\_modify, output, run start/stop, minimize info:**

The instantaneous values of the extended variables are written to
:doc:`binary restart files <restart>`.  Because the state of the random
number generator is not saved in restart files, this means you cannot
do "exact" restarts with this fix, where the simulation continues on
the same as if no restart had taken place. However, in a statistical
sense, a restarted simulation should produce the same behavior.
Note however that you should use a different seed each time you
restart, otherwise the same sequence of random numbers will be used
each time, which might lead to stochastic synchronization and
subtle artifacts in the sampling.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Langevin thermostatting to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".

Restrictions
""""""""""""


In order to perform constant-pressure simulations please use
:doc:`fix press/berendsen <fix_press_berendsen>`, rather than
:doc:`fix npt <fix_nh>`, to avoid duplicate integration of the
equations of motion.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix viscous <fix_viscous>`, :doc:`fix nvt <fix_nh>`, :doc:`pair\_style dpd/tstat <pair_dpd>`, :doc:`fix gld <fix_gld>`, :doc:`fix gle <fix_gle>`


----------


.. _Hijazi:



.. _Bussi3:

**(Hijazi)** M. Hijazi, D. M. Wilkins, M. Ceriotti, J. Chem. Phys. 148, 184109 (2018)


**(Bussi)** G. Bussi, M. Parrinello, Phs. Rev. E 75, 056707 (2007)


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
