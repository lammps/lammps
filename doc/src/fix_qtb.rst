.. index:: fix qtb

fix qtb command
===============

Syntax
""""""


.. parsed-literal::

   fix ID group-ID qtb keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* qtb = style name of this fix
* zero or more keyword/value pairs may be appended
* keyword = *temp* or *damp* or *seed* or *f\_max* or *N\_f*
  
  .. parsed-literal::
  
       *temp* value = target quantum temperature (temperature units)
       *damp* value = damping parameter (time units) inverse of friction *gamma*
       *seed* value = random number seed (positive integer)
       *f_max* value = upper cutoff frequency of the vibration spectrum (1/time units)
       *N_f* value = number of frequency bins (positive integer)



Examples
""""""""


.. code-block:: LAMMPS

   # (liquid methane modeled with the REAX force field, real units)
   fix 1 all nve
   fix 1 all qtb temp 110 damp 200 seed 35082 f_max 0.3 N_f 100
   # (quartz modeled with the BKS force field, metal units)
   fix 2 all nph iso 1.01325 1.01325 1
   fix 2 all qtb temp 300 damp 1 seed 47508 f_max 120.0 N_f 100

Description
"""""""""""

This command performs the quantum thermal bath scheme proposed by
:ref:`(Dammak) <Dammak>` to include self-consistent quantum nuclear effects,
when used in conjunction with the :doc:`fix nve <fix_nve>` or :doc:`fix nph <fix_nh>` commands.

Classical molecular dynamics simulation does not include any quantum
nuclear effect. Quantum treatment of the vibrational modes will
introduce zero point energy into the system, alter the energy power
spectrum and bias the heat capacity from the classical limit. Missing
all the quantum nuclear effects, classical MD cannot model systems at
temperatures lower than their classical limits. This effect is
especially important for materials with a large population of hydrogen
atoms and thus higher classical limits.

The equation of motion implemented by this command follows a Langevin
form:

.. math::

   m_i a_i = f_i + R_i - m_i\gamma v_i

Here :math:`m_i, a_i, f_i, R_i, \gamma, \textrm{and} v_i`
represent in this order mass, acceleration, force exerted by all other atoms, random
force, frictional coefficient (the inverse of damping parameter damp),
and velocity. The random force :math:`R_i` is "colored" so
that any vibrational mode with frequency :math:`\omega` will have a
temperature-sensitive energy :math:`\theta(\omega,T)` which
resembles the energy expectation for a quantum harmonic oscillator
with the same natural frequency:


.. math::

   \theta(\omega T) = \frac{\hbar}{2} + \hbar\omega \left[\exp(\frac{\hbar\omega}{k_B T})-1 \right]^{-1}

To efficiently generate the random forces, we employ the method
of :ref:`(Barrat) <Barrat>`, that circumvents the need to generate all
random forces for all times before the simulation. The memory
requirement of this approach is less demanding and independent
of the simulation duration. Since the total random force :math:`R_{tot}`
does not necessarily vanish for a finite number of atoms,
:math:`R_i` is replaced by :math:`R_i - \frac{R_{tot}}{N_{tot}}`
to avoid collective motion of the system.

The *temp* parameter sets the target quantum temperature. LAMMPS will
still have an output temperature in its thermo style. That is the
instantaneous classical temperature :math:`T^{cl}` derived from
the atom velocities at thermal equilibrium. A non-zero
:math:`T^{cl}` will be present even when the quantum
temperature approaches zero. This is associated with zero-point energy
at low temperatures.

.. math::

   T^{cl} = \sum \frac{m_i v_i^2}{3 N k_B}

The *damp* parameter is specified in time units, and it equals the
inverse of the frictional coefficient :math:`\gamma`. :math:`\gamma`
should be as small as possible but slightly larger than the timescale
of anharmonic coupling in the system which is about 10 ps to 100
ps. When :math:`\gamma` is too large, it gives an energy spectrum that
differs from the desired Bose-Einstein spectrum. When :math:`\gamma`
is too small, the quantum thermal bath coupling to the system will be
less significant than anharmonic effects, reducing to a classical
limit. We find that setting :math:`\gamma` between 5 THz and 1 THz
could be appropriate depending on the system.

The random number *seed* is a positive integer used to initiate a
Marsaglia random number generator. Each processor uses the input seed
to generate its own unique seed and its own stream of random
numbers. Thus the dynamics of the system will not be identical on two
runs on different numbers of processors.

The *f\_max* parameter truncate the noise frequency domain so that
vibrational modes with frequencies higher than *f\_max* will not be
modulated. If we denote :math:`\Delta t` as the time interval for the
MD integration, *f\_max* is always reset by the code to make
:math:`\alpha = (int)(2` *f\_max* :math:`\Delta t)^{-1}` a
positive integer and print out relative information. An appropriate
value for the cutoff frequency *f\_max* would be around 2~3 :math:`f_D`,
where :math:`f_D` is the Debye frequency.

The *N\_f* parameter is the frequency grid size, the number of points
from 0 to *f\_max* in the frequency domain that will be
sampled. 3*2\ *N\_f* per-atom random numbers are required
in the random force generation and there could be as many atoms as in
the whole simulation that can migrate into every individual
processor. A larger *N\_f* provides a more accurate sampling of the
spectrum while consumes more memory.  With fixed *f\_max* and
:math:`\gamma`, *N\_f* should be big enough to converge the classical
temperature :math:`T^{cl}` as a function of target quantum bath
temperature. Memory usage per processor could be from 10 to 100
Mbytes.

.. note::

   Unlike the :doc:`fix nvt <fix_nh>` command which performs
   Nose/Hoover thermostatting AND time integration, this fix does NOT
   perform time integration. It only modifies forces to a colored
   thermostat. Thus you must use a separate time integration fix, like
   :doc:`fix nve <fix_nve>` or :doc:`fix nph <fix_nh>` to actually
   update the velocities and positions of atoms (as shown in the
   examples). Likewise, this fix should not normally be used with
   other fixes or commands that also specify system temperatures ,
   e.g. :doc:`fix nvt <fix_nh>` and :doc:`fix temp/rescale
   <fix_temp_rescale>`.


----------


**Restart, fix\_modify, output, run start/stop, minimizie info:**

No information about this fix is written to :doc:`binary restart files
<restart>`.  Because the state of the random number generator is not
saved in restart files, this means you cannot do "exact" restarts with
this fix. However, in a statistical sense, a restarted simulation
should produce similar behaviors of the system.

This fix is not invoked during :doc:`energy minimization <minimize>`.


----------


Restrictions
""""""""""""


This fix style is part of the USER-QTB package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` doc page for more info.


----------


Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nph <fix_nh>`,
:doc:`fix langevin <fix_langevin>`, :doc:`fix qbmsst <fix_qbmsst>`


----------


Default
"""""""

The keyword defaults are temp = 300, damp = 1, seed = 880302,
f\_max=200.0 and N\_f = 100.


----------


.. _Dammak:



**(Dammak)** Dammak, Chalopin, Laroche, Hayoun, and Greffet, Phys Rev
Lett, 103, 190601 (2009).

.. _Barrat:



**(Barrat)** Barrat and Rodney, J. Stat. Phys, 144, 679 (2011).
