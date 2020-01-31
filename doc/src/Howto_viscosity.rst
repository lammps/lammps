Calculate viscosity
===================

The shear viscosity eta of a fluid can be measured in at least 5 ways
using various options in LAMMPS.  See the examples/VISCOSITY directory
for scripts that implement the 5 methods discussed here for a simple
Lennard-Jones fluid model.  Also, see the :doc:`Howto kappa <Howto_kappa>` doc page for an analogous discussion for
thermal conductivity.

Eta is a measure of the propensity of a fluid to transmit momentum in
a direction perpendicular to the direction of velocity or momentum
flow.  Alternatively it is the resistance the fluid has to being
sheared.  It is given by

J = -eta grad(Vstream)

where J is the momentum flux in units of momentum per area per time.
and grad(Vstream) is the spatial gradient of the velocity of the fluid
moving in another direction, normal to the area through which the
momentum flows.  Viscosity thus has units of pressure-time.

The first method is to perform a non-equilibrium MD (NEMD) simulation
by shearing the simulation box via the :doc:`fix deform <fix_deform>`
command, and using the :doc:`fix nvt/sllod <fix_nvt_sllod>` command to
thermostat the fluid via the SLLOD equations of motion.
Alternatively, as a second method, one or more moving walls can be
used to shear the fluid in between them, again with some kind of
thermostat that modifies only the thermal (non-shearing) components of
velocity to prevent the fluid from heating up.

.. note::

   A recent (2017) book by :ref:`(Daivis and Todd) <Daivis-viscosity>`
   discusses use of the SLLOD method and non-equilibrium MD (NEMD)
   thermostatting generally, for both simple and complex fluids,
   e.g. molecular systems.  The latter can be tricky to do correctly.

In both cases, the velocity profile setup in the fluid by this
procedure can be monitored by the :doc:`fix ave/chunk <fix_ave_chunk>`
command, which determines grad(Vstream) in the equation above.
E.g. the derivative in the y-direction of the Vx component of fluid
motion or grad(Vstream) = dVx/dy.  The Pxy off-diagonal component of
the pressure or stress tensor, as calculated by the :doc:`compute pressure <compute_pressure>` command, can also be monitored, which
is the J term in the equation above.  See the :doc:`Howto nemd <Howto_nemd>` doc page for details on NEMD simulations.

The third method is to perform a reverse non-equilibrium MD simulation
using the :doc:`fix viscosity <fix_viscosity>` command which implements
the rNEMD algorithm of Muller-Plathe.  Momentum in one dimension is
swapped between atoms in two different layers of the simulation box in
a different dimension.  This induces a velocity gradient which can be
monitored with the :doc:`fix ave/chunk <fix_ave_chunk>` command.
The fix tallies the cumulative momentum transfer that it performs.
See the :doc:`fix viscosity <fix_viscosity>` command for details.

The fourth method is based on the Green-Kubo (GK) formula which
relates the ensemble average of the auto-correlation of the
stress/pressure tensor to eta.  This can be done in a fully
equilibrated simulation which is in contrast to the two preceding
non-equilibrium methods, where momentum flows continuously through the
simulation box.

Here is an example input script that calculates the viscosity of
liquid Ar via the GK formalism:


.. parsed-literal::

   # Sample LAMMPS input script for viscosity of liquid Ar

   units       real
   variable    T equal 86.4956
   variable    V equal vol
   variable    dt equal 4.0
   variable    p equal 400     # correlation length
   variable    s equal 5       # sample interval
   variable    d equal $p\*$s   # dump interval

   # convert from LAMMPS real units to SI

   variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
   variable    atm2Pa equal 101325.0
   variable    A2m equal 1.0e-10
   variable    fs2s equal 1.0e-15
   variable    convert equal ${atm2Pa}\*${atm2Pa}\*${fs2s}\*${A2m}\*${A2m}\*${A2m}

   # setup problem

   dimension    3
   boundary     p p p
   lattice      fcc 5.376 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
   region       box block 0 4 0 4 0 4
   create_box   1 box
   create_atoms 1 box
   mass         1 39.948
   pair_style   lj/cut 13.0
   pair_coeff   \* \* 0.2381 3.405
   timestep     ${dt}
   thermo       $d

   # equilibration and thermalization

   velocity     all create $T 102486 mom yes rot yes dist gaussian
   fix          NVT all nvt temp $T $T 10 drag 0.2
   run          8000

   # viscosity calculation, switch to NVE if desired

   #unfix       NVT
   #fix         NVE all nve

   reset_timestep 0
   variable     pxy equal pxy
   variable     pxz equal pxz
   variable     pyz equal pyz
   fix          SS all ave/correlate $s $p $d &
                v_pxy v_pxz v_pyz type auto file S0St.dat ave running
   variable     scale equal ${convert}/(${kB}\*$T)\*$V\*$s\*${dt}
   variable     v11 equal trap(f_SS[3])\*${scale}
   variable     v22 equal trap(f_SS[4])\*${scale}
   variable     v33 equal trap(f_SS[5])\*${scale}
   thermo_style custom step temp press v_pxy v_pxz v_pyz v_v11 v_v22 v_v33
   run          100000
   variable     v equal (v_v11+v_v22+v_v33)/3.0
   variable     ndens equal count(all)/vol
   print        "average viscosity: $v [Pa.s] @ $T K, ${ndens} /A\^3"

The fifth method is related to the above Green-Kubo method,
but uses the Einstein formulation, analogous to the Einstein
mean-square-displacement formulation for self-diffusivity. The
time-integrated momentum fluxes play the role of Cartesian
coordinates, whose mean-square displacement increases linearly
with time at sufficiently long times.


----------


.. _Daivis-viscosity:



**(Daivis and Todd)** Daivis and Todd, Nonequilibrium Molecular Dynamics (book),
Cambridge University Press, https://doi.org/10.1017/9781139017848, (2017).


