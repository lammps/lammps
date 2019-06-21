.. index:: fix langevin/spin

fix langevin/spin command
=========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID langevin/spin T Tdamp seed

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/spin = style name of this fix command
* T = desired temperature of the bath (temperature units, K in metal units)
* Tdamp = transverse magnetic damping parameter (adim)
* seed = random number seed to use for white noise (positive integer)


Examples
""""""""


.. parsed-literal::

   fix 2 all langevin/spin 300.0 0.01 21

Description
"""""""""""

Apply a Langevin thermostat as described in :ref:`(Mayergoyz) <Mayergoyz1>` to the
magnetic spins associated to the atoms.
Used with :doc:`fix nve/spin <fix_nve_spin>`, this command performs
Brownian dynamics (BD).
A random torque and a transverse dissipation are applied to each spin i according to
the following stochastic differential equation:

.. math::

   :align: center

with lambda the transverse damping, and eta a random vector.
This equation is referred to as the stochastic Landau-Lifshitz-Gilbert (sLLG)
equation.

The components of eta are drawn from a Gaussian probability law. Their amplitude
is defined as a proportion of the temperature of the external thermostat T (in K
in metal units).

More details about this implementation are reported in :ref:`(Tranchida) <Tranchida2>`.

Note: due to the form of the sLLG equation, this fix has to be defined just
before the nve/spin fix (and after all other magnetic fixes).
As an example:


.. parsed-literal::

   fix 1 all precession/spin zeeman 0.01 0.0 0.0 1.0
   fix 2 all langevin/spin 300.0 0.01 21
   fix 3 all nve/spin lattice yes

is correct, but defining a force/spin command after the langevin/spin command
would give an error message.

Note: The random # *seed* must be a positive integer.  A Marsaglia random
number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  Because the state of the random number generator
is not saved in restart files, this means you cannot do "exact"
restarts with this fix, where the simulation continues on the same as
if no restart had taken place.  However, in a statistical sense, a
restarted simulation should produce the same behavior.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


The *langevin/spin* fix is part of the SPIN package.  This style is
only enabled if LAMMPS was built with this package.  See the :doc:`Build package <Build_package>` doc page for more info.

The numerical integration has to be performed with *fix nve/spin*
when *fix langevin/spin* is enabled.

This fix has to be the last defined magnetic fix before the time
integration fix (e.g. *fix nve/spin*\ ).

Related commands
""""""""""""""""

:doc:`fix nve/spin <fix_nve_spin>`, :doc:`fix precession/spin <fix_precession_spin>`

**Default:** none


----------


.. _Mayergoyz1:



**(Mayergoyz)** I.D. Mayergoyz, G. Bertotti, C. Serpico (2009). Elsevier (2009)

.. _Tranchida2:



**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
