.. index:: fix gjf

fix gjf command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID gjf Tstart Tstop damp seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* gjf = style name of this fix command
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* Tstart can be a variable (see below)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *vel* or *method*

  .. parsed-literal::

       *vel* value = *vfull* or *vhalf*
         *vfull* = use on-site velocity
         *vhalf* = use half-step velocity
       *method* value = *1-8*
         *1-8* = choose one of the many GJ formulations
         *7*   = requires input of additional scalar between 0 and 1

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 boundary gjf 10.0 10.0 1.0 699483
   fix 1 all gjf 10.0 100.0 100.0 48279 vel vfull method 4
   fix 2 all gjf 10.0 10.0 1.0 26488 method 7 0.95

Description
"""""""""""
.. versionadded:: 12Jun2025

Apply a Langevin thermostat as described in :ref:`(Gronbech-Jensen-2020) <Gronbech-Jensen-2020>`
to a group of atoms which models an interaction with a background
implicit solvent.  As described in the papers cited below, the GJ methods
provide exact diffusion, drift, and Boltzmann sampling for linear systems for
any time step within the stability limit. The purpose of this set of methods
is therefore to significantly improve statistical accuracy at longer time steps
compared to other thermostats.

The current implementation provides the user with the option to output
the velocity in one of two forms: *vfull* or *vhalf*. The option *vhalf*
outputs the 2GJ half-step velocity given in :ref:`Gronbech Jensen/Gronbech-Jensen
<Gronbech-Jensen-2019>`; for linear systems, this velocity is shown to not
have any statistical errors for any stable time step. The option *vfull*
outputs the on-site velocity given in :ref:`Gronbech-Jensen/Farago
<Gronbech-Jensen-Farago>`; this velocity is shown to be systematically lower
than the target temperature by a small amount, which grows
quadratically with the timestep. An overview of statistically correct Boltzmann
and Maxwell-Boltzmann sampling of true on-site and true half-step velocities is
given in :ref:`Gronbech-Jensen-2020 <Gronbech-Jensen-2020>`.

This fix allows the use of several GJ methods as listed in :ref:`Gronbech-Jensen-2020 <Gronbech-Jensen-2020>`.
The GJ-VII method is described in :ref:`Finkelstein <Finkelstein>` and GJ-VIII
is described in :ref:`Gronbech-Jensen-2024 <Gronbech-Jensen-2024>`.
The implementation follows the splitting form provided in Eqs. (24) and (25)
in :ref:`Gronbech-Jensen-2024 <Gronbech-Jensen-2024>`, including the application
of Gaussian noise values, per the description in
:ref:`Gronbech-Jensen-2023 <Gronbech-Jensen-2023>`.


.. note::

   Unlike the :doc:`fix langevin <fix_langevin>` command which performs force
   modifications only, this fix performs thermostatting and time integration.
   Thus you no longer need a separate time integration fix, like :doc:`fix nve <fix_nve>`.

See the :doc:`Howto thermostat <Howto_thermostat>` page for
a discussion of different ways to compute temperature and perform
thermostatting.

The desired temperature at each timestep is a ramped value during the
run from *Tstart* to *Tstop*\ .

*Tstart* can be specified as an equal-style or atom-style
:doc:`variable <variable>`.  In this case, the *Tstop* setting is
ignored.  If the value is a variable, it should be specified as
v_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
target temperature.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent temperature.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent
temperature with optional time-dependence as well.

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
a bias.

The *damp* parameter is specified in time units and determines how
rapidly the temperature is relaxed.  For example, a value of 100.0 means
to relax the temperature in a timespan of (roughly) 100 time units
(:math:`\tau` or fs or ps - see the :doc:`units <units>` command).  The
damp factor can be thought of as inversely related to the viscosity of
the solvent.  I.e. a small relaxation time implies a high-viscosity
solvent and vice versa.  See the discussion about :math:`\gamma` and
viscosity in the documentation for the :doc:`fix viscous <fix_viscous>`
command for more details.

The random # *seed* must be a positive integer.  A Marsaglia random
number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.

----------

The keyword/value option pairs are used in the following ways.

The keyword *vel* determines which velocity is used to determine
quantities of interest in the simulation.

The keyword *method* selects one of the eight GJ-methods implemented in LAMMPS.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
Because the state of the random number generator is not saved in restart files,
this means you cannot do "exact" restarts with this fix, where the simulation
continues on the same as if no restart had taken place.  However, in a
statistical sense, a restarted simulation should produce the same behavior.
Additionally, the GJ methods implement noise exclusively within each time step
(unlike the BBK thermostat of the fix-langevin). The restart is done with
either vfull or vhalf velocity output for as long as the choice of vfull/vhalf
is the same for the simulation as it is in the restart file.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is not compatible with run_style respa. It is not compatible with
accelerated packages such as KOKKOS.

Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`, :doc:`fix nvt <fix_nh>`

Default
"""""""

The option defaults are vel = vhalf, method = 1.

----------

.. _Gronbech-Jensen-2020:

**(Gronbech-Jensen-2020)** Gronbech-Jensen, Mol Phys 118, e1662506 (2020).

.. _Gronbech-Jensen-2019:

**(Gronbech Jensen/Gronbech-Jensen)** Gronbech Jensen and Gronbech-Jensen, Mol Phys, 117, 2511 (2019)

.. _Gronbech-Jensen-Farago:

**(Gronbech-Jensen/Farago)** Gronbech-Jensen and Farago, Mol Phys, 111, 983 (2013).

.. _Finkelstein:

**(Finkelstein)** Finkelstein, Cheng, Florin, Seibold, Gronbech-Jensen, J. Chem. Phys., 155, 18 (2021)

.. _Gronbech-Jensen-2024:

**(Gronbech-Jensen-2024)** Gronbech-Jensen, J. Stat. Phys. 191, 137 (2024).

.. _Gronbech-Jensen-2023:

**(Gronbech-Jensen-2023)** Gronbech-Jensen, J. Stat. Phys. 190, 96 (2023).
