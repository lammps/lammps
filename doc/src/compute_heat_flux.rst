.. index:: compute heat/flux

compute heat/flux command
=========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID heat/flux ke-ID pe-ID stress-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* heat/flux = style name of this compute command
* ke-ID = ID of a compute that calculates per-atom kinetic energy
* pe-ID = ID of a compute that calculates per-atom potential energy
* stress-ID = ID of a compute that calculates per-atom stress

Examples
""""""""

.. code-block:: LAMMPS

   compute myFlux all heat/flux myKE myPE myStress

Description
"""""""""""

Define a computation that calculates the heat flux vector based on
contributions from atoms in the specified group.  This can be used by
itself to measure the heat flux through a set of atoms (e.g. a region
between two thermostatted reservoirs held at different temperatures),
or to calculate a thermal conductivity using the equilibrium
Green-Kubo formalism.

For other non-equilibrium ways to compute a thermal conductivity, see
the :doc:`Howto kappa <Howto_kappa>` doc page..  These include use of
the :doc:`fix thermal/conductivity <fix_thermal_conductivity>` command
for the Muller-Plathe method.  Or the :doc:`fix heat <fix_heat>` command
which can add or subtract heat from groups of atoms.

The compute takes three arguments which are IDs of other
:doc:`computes <compute>`.  One calculates per-atom kinetic energy
(\ *ke-ID*\ ), one calculates per-atom potential energy (\ *pe-ID)*, and the
third calculates per-atom stress (\ *stress-ID*\ ).

.. note::

   These other computes should provide values for all the atoms in
   the group this compute specifies.  That means the other computes could
   use the same group as this compute, or they can just use group "all"
   (or any group whose atoms are superset of the atoms in this compute's
   group).  LAMMPS does not check for this.

In case of two-body interactions, the heat flux is defined as:

.. math::
   \mathbf{J} &= \frac{1}{V} \left[ \sum_i e_i \mathbf{v}_i - \sum_{i} \mathbf{S}_{i} \mathbf{v}_i \right] \\
   &= \frac{1}{V} \left[ \sum_i e_i \mathbf{v}_i + \sum_{i<j} \left( \mathbf{F}_{ij} \cdot \mathbf{v}_j \right) \mathbf{r}_{ij} \right] \\
   &= \frac{1}{V} \left[ \sum_i e_i \mathbf{v}_i + \frac{1}{2} \sum_{i<j} \left( \mathbf{F}_{ij} \cdot \left(\mathbf{v}_i + \mathbf{v}_j \right)  \right) \mathbf{r}_{ij} \right]

:math:`e_i` in the first term of the equation
is the per-atom energy (potential and kinetic).
This is calculated by the computes *ke-ID*
and *pe-ID*. :math:`\mathbf{S}_i` in the second term is the
per-atom stress tensor calculated by the compute *stress-ID*.
See :doc:`compute stress/atom <compute_stress_atom>`
and :doc:`compute centroid/stress/atom <compute_stress_atom>`
for possible definitions of atomic stress :math:`\mathbf{S}_i`
in the case of bonded and many-body interactions.
The tensor multiplies :math:`\mathbf{v}_i` as a 3x3 matrix-vector multiply
to yield a vector.
Note that as discussed below, the 1/:math:`{V}` scaling factor in the
equation for :math:`\mathbf{J}` is NOT included in the calculation performed by
these computes; you need to add it for a volume appropriate to the atoms
included in the calculation.

.. note::

   The :doc:`compute pe/atom <compute_pe_atom>` and
   :doc:`compute stress/atom <compute_stress_atom>`
   commands have options for which
   terms to include in their calculation (pair, bond, etc).  The heat
   flux calculation will thus include exactly the same terms. Normally
   you should use :doc:`compute stress/atom virial <compute_stress_atom>`
   or :doc:`compute centroid/stress/atom virial <compute_stress_atom>`
   so as not to include a kinetic energy term in the heat flux.

.. warning::

   The compute *heat/flux* has been reported to produce unphysical
   values for angle, dihedral and improper contributions
   when used with :doc:`compute stress/atom <compute_stress_atom>`,
   as discussed in :ref:`(Surblys) <Surblys2>` and :ref:`(Boone) <Boone>`.
   You are strongly advised to
   use :doc:`compute centroid/stress/atom <compute_stress_atom>`,
   which has been implemented specifically for such cases.

The Green-Kubo formulas relate the ensemble average of the
auto-correlation of the heat flux :math:`\mathbf{J}`
to the thermal conductivity :math:`\kappa`:

.. math::
   \kappa  = \frac{V}{k_B T^2} \int_0^\infty \langle J_x(0)  J_x(t) \rangle \, \mathrm{d} t = \frac{V}{3 k_B T^2} \int_0^\infty \langle \mathbf{J}(0) \cdot  \mathbf{J}(t)  \rangle \, \mathrm{d}t

----------

The heat flux can be output every so many timesteps (e.g. via the
:doc:`thermo_style custom <thermo_style>` command).  Then as a
post-processing operation, an auto-correlation can be performed, its
integral estimated, and the Green-Kubo formula above evaluated.

The :doc:`fix ave/correlate <fix_ave_correlate>` command can calculate
the auto-correlation.  The trap() function in the
:doc:`variable <variable>` command can calculate the integral.

An example LAMMPS input script for solid Ar is appended below.  The
result should be: average conductivity ~0.29 in W/mK.

----------

Output info
"""""""""""

This compute calculates a global vector of length 6.
The first 3 components are the :math:`x`, :math:`y`, :math:`z`
components of the full heat flux vector,
i.e. (:math:`J_x`, :math:`J_y`, :math:`J_z`).
The next 3 components are the :math:`x`, :math:`y`, :math:`z` components
of just the convective portion of the flux, i.e. the
first term in the equation for :math:`\mathbf{J}`.
Each component can be
accessed by indices 1-6. These values can be used by any command that
uses global vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The vector values calculated by this compute are "extensive", meaning
they scale with the number of atoms in the simulation.  They can be
divided by the appropriate volume to get a flux, which would then be
an "intensive" value, meaning independent of the number of atoms in
the simulation.  Note that if the compute is "all", then the
appropriate volume to divide by is the simulation box volume.
However, if a sub-group is used, it should be the volume containing
those atoms.

The vector values will be in energy\*velocity :doc:`units <units>`.  Once
divided by a volume the units will be that of flux, namely
energy/area/time :doc:`units <units>`

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix thermal/conductivity <fix_thermal_conductivity>`,
:doc:`fix ave/correlate <fix_ave_correlate>`,
:doc:`variable <variable>`

Default
"""""""

none

----------

.. code-block:: LAMMPS

   # Sample LAMMPS input script for thermal conductivity of solid Ar

   units       real
   variable    T equal 70
   variable    V equal vol
   variable    dt equal 4.0
   variable    p equal 200     # correlation length
   variable    s equal 10      # sample interval
   variable    d equal $p*$s   # dump interval

   # convert from LAMMPS real units to SI

   variable    kB equal 1.3806504e-23    # [J/K] Boltzmann
   variable    kCal2J equal 4186.0/6.02214e23
   variable    A2m equal 1.0e-10
   variable    fs2s equal 1.0e-15
   variable    convert equal ${kCal2J}*${kCal2J}/${fs2s}/${A2m}

   # setup problem

   dimension    3
   boundary     p p p
   lattice      fcc 5.376 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
   region       box block 0 4 0 4 0 4
   create_box   1 box
   create_atoms 1 box
   mass         1 39.948
   pair_style   lj/cut 13.0
   pair_coeff   * * 0.2381 3.405
   timestep     ${dt}
   thermo       $d

   # equilibration and thermalization

   velocity     all create $T 102486 mom yes rot yes dist gaussian
   fix          NVT all nvt temp $T $T 10 drag 0.2
   run          8000

   # thermal conductivity calculation, switch to NVE if desired

   #unfix       NVT
   #fix         NVE all nve

   reset_timestep 0
   compute      myKE all ke/atom
   compute      myPE all pe/atom
   compute      myStress all stress/atom NULL virial
   compute      flux all heat/flux myKE myPE myStress
   variable     Jx equal c_flux[1]/vol
   variable     Jy equal c_flux[2]/vol
   variable     Jz equal c_flux[3]/vol
   fix          JJ all ave/correlate $s $p $d &
                c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running
   variable     scale equal ${convert}/${kB}/$T/$T/$V*$s*${dt}
   variable     k11 equal trap(f_JJ[3])*${scale}
   variable     k22 equal trap(f_JJ[4])*${scale}
   variable     k33 equal trap(f_JJ[5])*${scale}
   thermo_style custom step temp v_Jx v_Jy v_Jz v_k11 v_k22 v_k33
   run          100000
   variable     k equal (v_k11+v_k22+v_k33)/3.0
   variable     ndens equal count(all)/vol
   print        "average conductivity: $k[W/mK] @ $T K, ${ndens} /A\^3"

----------

.. _Surblys2:

**(Surblys)** Surblys, Matsubara, Kikugawa, Ohara, Phys Rev E, 99, 051301(R) (2019).

.. _Boone:

**(Boone)** Boone, Babaei, Wilmer, J Chem Theory Comput, 15, 5579--5587 (2019).
