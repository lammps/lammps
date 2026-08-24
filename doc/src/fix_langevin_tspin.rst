.. index:: fix langevin/tspin

fix langevin/tspin command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID langevin/tspin Tstart Tstop Tdamp seed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/tspin = style name of this fix command
* Tstart, Tstop = desired spin temperature at start and end of the run
  (temperature units, K in metal units)
* Tdamp = spin temperature damping parameter (time units)
* seed = random number seed to use for the white noise (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *zero*

  .. parsed-literal::

       *zero* value = *no* or *yes*
         no = the net random force on the spin system need not be zero
         yes = subtract the mean random force so that it sums to zero

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all langevin/tspin 300.0 300.0 0.05 48279
   fix 2 all langevin/tspin 10.0 600.0 0.1 48279 zero yes

Description
"""""""""""

.. versionadded:: TBD

Apply a Langevin thermostat to the spin velocities of inertial spin
dynamics.  A friction term and a random force are added to the magnetic
force of every magnetic atom in the group,

.. math::

   \vec{F}^{m}_{i} \rightarrow \vec{F}^{m}_{i}
   - \frac{m_s}{T_{damp}} \vec{v}^{s}_i
   + \sqrt{\frac{2 m_s k_B T_s}{T_{damp}\, dt}}\, \vec{R}_i

where :math:`\vec{R}_i` is a vector of uniform random numbers with zero
mean and unit variance.  The two terms together make the spin
velocities sample a Maxwell-Boltzmann distribution at the target spin
temperature, which is ramped linearly from Tstart to Tstop over the
course of the run.

This fix does not perform time integration; combine it with :doc:`fix
nve/tspin <fix_nve_tspin>`.  The warning in that page about the need
for a longitudinal potential on the spin modulus applies here as well.

This thermostat acts on the spin momenta and is unrelated to
:doc:`fix langevin/spin <fix_langevin_spin>`, which adds transverse
Gilbert damping and a fluctuating field to the fixed-modulus
Landau-Lifshitz-Gilbert equation.  The two fixes should not be combined.

Note: the random number *seed* must be a positive integer.  A Marsaglia
random number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  Because the state of the random number generator is not
saved in restart files, this means you cannot do "exact" restarts with
this fix, where the simulation continues on the same as if no restart
had taken place.  However, in a statistical sense, a restarted
simulation should produce the same behavior.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The *langevin/tspin* fix is part of the SPIN package.  This style is
only enabled if LAMMPS was built with this package.  See the
:doc:`Build package <Build_package>` page for more info.

This fix requires :doc:`atom_style tspin <atom_style>`.  The time
integration has to be performed by :doc:`fix nve/tspin
<fix_nve_tspin>`.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`fix langevin <fix_langevin>`

Default
""""""""

The default value for the *zero* keyword is *no*.
