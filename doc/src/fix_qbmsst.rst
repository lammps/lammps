.. index:: fix qbmsst

fix qbmsst command
==================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID qbmsst dir shockvel keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* qbmsst = style name of this fix
* dir = *x* or *y* or *z*
* shockvel = shock velocity (strictly positive, velocity units)
* zero or more keyword/value pairs may be appended
* keyword = *q* or *mu* or *p0* or *v0* or *e0* or *tscale* or *damp* or *seed*\ or *f\_max* or *N\_f* or *eta* or *beta* or *T\_init*
  
  .. parsed-literal::
  
       *q* value = cell mass-like parameter (mass\^2/distance\^4 units)
       *mu* value = artificial viscosity (mass/distance/time units)
       *p0* value = initial pressure in the shock equations (pressure units)
       *v0* value = initial simulation cell volume in the shock equations (distance\^3 units)
       *e0* value = initial total energy (energy units)
       *tscale* value = reduction in initial temperature (unitless fraction between 0.0 and 1.0)
       *damp* value = damping parameter (time units) inverse of friction <i>&gamma;</i>
       *seed* value = random number seed (positive integer)
       *f_max* value = upper cutoff frequency of the vibration spectrum (1/time units)
       *N_f* value = number of frequency bins (positive integer)
       *eta* value = coupling constant between the shock system and the quantum thermal bath (positive unitless)
       *beta* value = the quantum temperature is updated every beta time steps (positive integer)
       *T_init* value = quantum temperature for the initial state (temperature units)



Examples
""""""""


.. parsed-literal::

   fix 1 all qbmsst z 0.122 q 25 mu 0.9 tscale 0.01 damp 200 seed 35082 f_max 0.3 N_f 100 eta 1 beta 400 T_init 110 (liquid methane modeled with the REAX force field, real units)
   fix 2 all qbmsst z 72 q 40 tscale 0.05 damp 1 seed 47508 f_max 120.0 N_f 100 eta 1.0 beta 500 T_init 300 (quartz modeled with the BKS force field, metal units)

Two example input scripts are given, including shocked
:math:`\alpha-\mathrm{quartz}` and shocked liquid methane.
The input script first equilibrate an initial state with the quantum
thermal bath at the target temperature and then apply the qbmsst to
simulate shock compression with quantum nuclear correction.  The
following two figures plot related quantities for shocked alpha quartz.

.. image:: JPG/qbmsst_init.jpg
   :align: center

Figure 1. Classical temperature <i>T</i><sup>cl</sup> = &sum;
<i>m<sub>i</sub>v<sub>i</sub><sup>2</sup>/3Nk</i><sub>B</sub> vs. time
for coupling the alpha quartz initial state with the quantum thermal
bath at target quantum temperature <i>T</i><sup>qm</sup> = 300 K. The
NpH ensemble is used for time integration while QTB provides the
colored random force. <i>T</i><sup>cl</sup> converges at the timescale
of *damp* which is set to be 1 ps.

.. image:: JPG/qbmsst_shock.jpg
   :align: center

Figure 2. Quantum temperature and pressure vs. time for simulating
shocked alpha quartz with the QBMSST. The shock propagates along the z
direction. Restart of the QBMSST command is demonstrated in the
example input script. Thermodynamic quantities stay continuous before
and after the restart.

Description
"""""""""""

This command performs the Quantum-Bath coupled Multi-Scale Shock
Technique (QBMSST) integration. See :ref:`(Qi) <Qi>` for a detailed
description of this method.  The QBMSST provides description of the
thermodynamics and kinetics of shock processes while incorporating
quantum nuclear effects.  The *shockvel* setting determines the steady
shock velocity that will be simulated along direction *dir*\ .

Quantum nuclear effects :doc:`(fix qtb) <fix_qtb>` can be crucial
especially when the temperature of the initial state is below the
classical limit or there is a great change in the zero point energies
between the initial and final states. Theoretical post processing
quantum corrections of shock compressed water and methane have been
reported as much as 30% of the temperatures :ref:`(Goldman) <Goldman1>`.  A
self-consistent method that couples the shock to a quantum thermal
bath described by a colored noise Langevin thermostat has been
developed by Qi et al :ref:`(Qi) <Qi>` and applied to shocked methane.  The
onset of chemistry is reported to be at a pressure on the shock
Hugoniot that is 40% lower than observed with classical molecular
dynamics.

It is highly recommended that the system be already in an equilibrium
state with a quantum thermal bath at temperature of *T\_init*.  The fix
command :doc:`fix qtb <fix_qtb>` at constant temperature *T\_init* could
be used before applying this command to introduce self-consistent
quantum nuclear effects into the initial state.

The parameters *q*\ , *mu*\ , *e0*\ , *p0*\ , *v0* and *tscale* are described
in the command :doc:`fix msst <fix_msst>`. The values of *e0*\ , *p0*\ , or
*v0* will be calculated on the first step if not specified.  The
parameter of *damp*\ , *f\_max*, and *N\_f* are described in the command
:doc:`fix qtb <fix_qtb>`.

The fix qbmsst command couples the shock system to a quantum thermal
bath with a rate that is proportional to the change of the total
energy of the shock system, <i>etot</i> - <i>etot</i><sub>0</sub>.
Here <i>etot</i> consists of both the system energy and a thermal
term, see :ref:`(Qi) <Qi>`, and <i>etot</i><sub>0</sub> = *e0* is the
initial total energy.

The *eta* (<i>&eta;</i>) parameter is a unitless coupling constant
between the shock system and the quantum thermal bath. A small *eta*
value cannot adjust the quantum temperature fast enough during the
temperature ramping period of shock compression while large *eta*
leads to big temperature oscillation. A value of *eta* between 0.3 and
1 is usually appropriate for simulating most systems under shock
compression. We observe that different values of *eta* lead to almost
the same final thermodynamic state behind the shock, as expected.

The quantum temperature is updated every *beta* (<i>&beta;</i>) steps
with an integration time interval *beta* times longer than the
simulation time step. In that case, <i>etot</i> is taken as its
average over the past *beta* steps. The temperature of the quantum
thermal bath <i>T</i><sup>qm</sup> changes dynamically according to
the following equation where &Delta;<i>t</i> is the MD time step and
<i>&gamma;</i> is the friction constant which is equal to the inverse
of the *damp* parameter.

.. raw:: html

   <center><font size="4"> <i>dT</i><sup>qm</sup>/<i>dt =
   &gamma;&eta;</i>&sum;<i><sup>&beta;</sup><sub>l =
   1</sub></i>[<i>etot</i>(<i>t-l</i>&Delta;<i>t</i>)-<i>etot</i><sub>0</sub>]/<i>3&beta;Nk</i><sub>B</sub>
   </font></center>

The parameter *T\_init* is the initial temperature of the quantum
thermal bath and the system before shock loading.

For all pressure styles, the simulation box stays orthorhombic in
shape. Parrinello-Rahman boundary conditions (tilted box) are
supported by LAMMPS, but are not implemented for QBMSST.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

Because the state of the random number generator is not written to
:doc:`binary restart files <restart>`, this fix cannot be restarted
"exactly" in an uninterrupted fashion. However, in a statistical
sense, a restarted simulation should produce similar behaviors of the
system as if it is not interrupted.  To achieve such a restart, one
should write explicitly the same value for *q*\ , *mu*\ , *damp*\ , *f\_max*,
*N\_f*, *eta*\ , and *beta* and set *tscale* = 0 if the system is
compressed during the first run.

The progress of the QBMSST can be monitored by printing the global
scalar and global vector quantities computed by the fix.  The global
vector contains five values in this order:

[\ *dhugoniot*\ , *drayleigh*\ , *lagrangian\_speed*, *lagrangian\_position*,
*quantum\_temperature*]

1. *dhugoniot* is the departure from the Hugoniot (temperature units).
2. *drayleigh* is the departure from the Rayleigh line (pressure units).
3. *lagrangian\_speed* is the laboratory-frame Lagrangian speed (particle velocity) of the computational cell (velocity units).
4. *lagrangian\_position* is the computational cell position in the reference frame moving at the shock speed. This is the distance of the computational cell behind the shock front.
5. *quantum\_temperature* is the temperature of the quantum thermal bath <i>T</i><sup>qm</sup>.

To print these quantities to the log file with descriptive column
headers, the following LAMMPS commands are suggested. Here the
:doc:`fix_modify <fix_modify>` energy command is also enabled to allow
the thermo keyword *etotal* to print the quantity <i>etot</i>.  See
also the :doc:`thermo_style <thermo_style>` command.


.. parsed-literal::

   fix             fix_id all msst z
   fix_modify      fix_id energy yes
   variable        dhug    equal f_fix_id[1]
   variable        dray    equal f_fix_id[2]
   variable        lgr_vel equal f_fix_id[3]
   variable        lgr_pos equal f_fix_id[4]
   variable        T_qm    equal f_fix_id[5]
   thermo_style    custom  step temp ke pe lz pzz etotal v_dhug v_dray v_lgr_vel v_lgr_pos v_T_qm f_fix_id

The global scalar under the entry f\_fix\_id is the quantity of thermo
energy as an extra part of <i>etot</i>. This global scalar and the
vector of 5 quantities can be accessed by various :doc:`output commands <Howto_output>`. It is worth noting that the temp keyword
under the :doc:`thermo_style <thermo_style>` command print the
instantaneous classical temperature <i>T</i><sup>cl</sup> as described
in the command :doc:`fix qtb <fix_qtb>`.


----------


Restrictions
""""""""""""


This fix style is part of the USER-QTB package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

All cell dimensions must be periodic. This fix can not be used with a
triclinic cell.  The QBMSST fix has been tested only for the group-ID
all.


----------


Related commands
""""""""""""""""

:doc:`fix qtb <fix_qtb>`, :doc:`fix msst <fix_msst>`


----------


Default
"""""""

The keyword defaults are q = 10, mu = 0, tscale = 0.01, damp = 1, seed
= 880302, f\_max = 200.0, N\_f = 100, eta = 1.0, beta = 100, and
T\_init=300.0. e0, p0, and v0 are calculated on the first step.


----------


.. _Goldman1:



**(Goldman)** Goldman, Reed and Fried, J. Chem. Phys. 131, 204103 (2009)

.. _Qi:



**(Qi)** Qi and Reed, J. Phys. Chem. A 116, 10451 (2012).
