.. index:: fix msst

fix msst command
================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID msst dir shockvel keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* msst = style name of this fix
* dir = *x* or *y* or *z*
* shockvel = shock velocity (strictly positive, distance/time units)
* zero or more keyword value pairs may be appended
* keyword = *q* or *mu* or *p0* or *v0* or *e0* or *tscale* or *beta* or *dftb*
  
  .. parsed-literal::
  
       *q* value = cell mass-like parameter (mass\^2/distance\^4 units)
       *mu* value = artificial viscosity (mass/length/time units)
       *p0* value = initial pressure in the shock equations (pressure units)
       *v0* value = initial simulation cell volume in the shock equations (distance\^3 units)
       *e0* value = initial total energy (energy units)
       *tscale* value = reduction in initial temperature (unitless fraction between 0.0 and 1.0)
       *dftb* value = *yes* or *no* for whether using MSST in conjunction with DFTB+
       *beta* value = scale factor for improved energy conservation



Examples
""""""""


.. parsed-literal::

   fix 1 all msst y 100.0 q 1.0e5 mu 1.0e5
   fix 2 all msst z 50.0 q 1.0e4 mu 1.0e4  v0 4.3419e+03 p0 3.7797e+03 e0 -9.72360e+02 tscale 0.01
   fix 1 all msst y 100.0 q 1.0e5 mu 1.0e5 dftb yes beta 0.5

Description
"""""""""""

This command performs the Multi-Scale Shock Technique (MSST)
integration to update positions and velocities each timestep to mimic
a compressive shock wave passing over the system. See :ref:`(Reed) <Reed>`
for a detailed description of this method.  The MSST varies the cell
volume and temperature in such a way as to restrain the system to the
shock Hugoniot and the Rayleigh line. These restraints correspond to
the macroscopic conservation laws dictated by a shock
front. *shockvel* determines the steady shock velocity that will be
simulated.

To perform a simulation, choose a value of *q* that provides volume
compression on the timescale of 100 fs to 1 ps.  If the volume is not
compressing, either the shock speed is chosen to be below the material
sound speed or *p0* has been chosen inaccurately.  Volume compression
at the start can be sped up by using a non-zero value of *tscale*\ . Use
the smallest value of *tscale* that results in compression.

Under some special high-symmetry conditions, the pressure (volume)
and/or temperature of the system may oscillate for many cycles even
with an appropriate choice of mass-like parameter *q*\ . Such
oscillations have physical significance in some cases.  The optional
*mu* keyword adds an artificial viscosity that helps break the system
symmetry to equilibrate to the shock Hugoniot and Rayleigh line more
rapidly in such cases.

The keyword *tscale* is a factor between 0 and 1 that determines what
fraction of thermal kinetic energy is converted to compressive strain
kinetic energy at the start of the simulation.  Setting this parameter
to a non-zero value may assist in compression at the start of
simulations where it is slow to occur.

If keywords *e0*\ , *p0*\ ,or *v0* are not supplied, these quantities will
be calculated on the first step, after the energy specified by
*tscale* is removed.  The value of *e0* is not used in the dynamical
equations, but is used in calculating the deviation from the Hugoniot.

The keyword *beta* is a scaling term that can be added to the MSST
ionic equations of motion to account for drift in the conserved
quantity during long timescale simulations, similar to a Berendsen
thermostat. See :ref:`(Reed) <Reed>` and :ref:`(Goldman) <Goldman2>` for more
details.  The value of *beta* must be between 0.0 and 1.0 inclusive.
A value of 0.0 means no contribution, a value of 1.0 means a full
contribution.

Values of shockvel less than a critical value determined by the
material response will not have compressive solutions. This will be
reflected in lack of significant change of the volume in the MSST.

For all pressure styles, the simulation box stays orthogonal in shape.
Parrinello-Rahman boundary conditions (tilted box) are supported by
LAMMPS, but are not implemented for MSST.

This fix computes a temperature and pressure and potential energy each
timestep. To do this, the fix creates its own computes of style "temp"
"pressure", and "pe", as if these commands had been issued:


.. parsed-literal::

   compute fix-ID_MSST_temp all temp
   compute fix-ID_MSST_press all pressure fix-ID_MSST_temp

   compute fix-ID_MSST_pe all pe

See the :doc:`compute temp <compute_temp>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + "_MSST\_temp`or <MSST_press">`_
or "_MSST\_pe".  The group for the new computes is "all".


----------


The *dftb* keyword is to allow this fix to be used when LAMMPS is
being driven by DFTB+, a density-functional tight-binding code. If the
keyword *dftb* is used with a value of *yes*\ , then the MSST equations
are altered to account for the electron entropy contribution to the
Hugonio relations and total energy.  See :ref:`(Reed2) <Reed2>` and
:ref:`(Goldman) <Goldman2>` for details on this contribution.  In this case,
you must define a :doc:`fix external <fix_external>` command in your
input script, which is used to callback to DFTB+ during the LAMMPS
timestepping.  DFTB+ will communicate its info to LAMMPS via that fix.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the state of all internal variables to :doc:`binary restart files <restart>`.  See the :doc:`read_restart <read_restart>` command
for info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The progress of the MSST can be monitored by printing the global
scalar and global vector quantities computed by the fix.

The scalar is the cumulative energy change due to the fix. This is
also the energy added to the potential energy by the
:doc:`fix_modify <fix_modify>` *energy* command.  With this command, the
thermo keyword *etotal* prints the conserved quantity of the MSST
dynamic equations. This can be used to test if the MD timestep is
sufficiently small for accurate integration of the dynamic
equations. See also :doc:`thermo_style <thermo_style>` command.

The global vector contains four values in this order:

[\ *dhugoniot*\ , *drayleigh*\ , *lagrangian\_speed*, *lagrangian\_position*]

1. *dhugoniot* is the departure from the Hugoniot (temperature units).
2. *drayleigh* is the departure from the Rayleigh line (pressure units).
3. *lagrangian\_speed* is the laboratory-frame Lagrangian speed (particle velocity) of the computational cell (velocity units).
4. *lagrangian\_position* is the computational cell position in the reference frame moving at the shock speed. This is usually a good estimate of distance of the computational cell behind the shock front.

To print these quantities to the log file with descriptive column
headers, the following LAMMPS commands are suggested:


.. parsed-literal::

   fix              msst all msst z
   fix_modify       msst energy yes
   variable dhug    equal f_msst[1]
   variable dray    equal f_msst[2]
   variable lgr_vel equal f_msst[3]
   variable lgr_pos equal f_msst[4]
   thermo_style     custom step temp ke pe lz pzz etotal v_dhug v_dray v_lgr_vel v_lgr_pos f_msst

These fixes compute a global scalar and a global vector of 4
quantities, which can be accessed by various :doc:`output commands <Howto_output>`.  The scalar values calculated by this fix
are "extensive"; the vector values are "intensive".

Restrictions
""""""""""""


This fix style is part of the SHOCK package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

All cell dimensions must be periodic. This fix can not be used with a
triclinic cell.  The MSST fix has been tested only for the group-ID
all.

Related commands
""""""""""""""""

:doc:`fix nphug <fix_nphug>`, :doc:`fix deform <fix_deform>`

Default
"""""""

The keyword defaults are q = 10, mu = 0, tscale = 0.01, dftb = no,
beta = 0.0.  Note that p0, v0, and e0 are calculated on the first
timestep.


----------


.. _Reed:



**(Reed)** Reed, Fried, and Joannopoulos, Phys. Rev. Lett., 90, 235503
(2003).

.. _Reed2:



**(Reed2)** Reed, J. Phys. Chem. C, 116, 2205 (2012).

.. _Goldman2:



**(Goldman)** Goldman, Srinivasan, Hamel, Fried, Gaus, and Elstner,
J. Phys. Chem. C, 117, 7885 (2013).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
