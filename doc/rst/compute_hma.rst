.. index:: compute hma

compute hma command
===================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID hma temp-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* hma = style name of this compute command
* temp-ID = ID of fix that specifies the set temperature during canonical simulation
* keyword = *anharmonic* *u* *p Pharm* *cv*

.. parsed-literal::

     *anharmonic* = compute will return anharmonic property values
     *u* = compute will return potential energy
     *p* = compute will return pressure.  the following keyword must be the difference between the harmonic pressure and lattice pressure as described below
     *cv* = compute will return the heat capacity



Examples
""""""""


.. parsed-literal::

   compute 2 all hma 1 u
   compute 2 all hma 1 anharmonic u p 0.9
   compute 2 all hma 1 u cv

Description
"""""""""""

Define a computation that calculates the properties of a solid (potential
energy, pressure or heat capacity), using the harmonically-mapped averaging
(HMA) method. 
This command yields much higher precision than the equivalent compute commands
(:doc:`compute pe <compute_pe>`, :doc:`compute pressure <compute_pressure>`, etc.)
commands during a canonical simulation of an atomic crystal. Specifically,
near melting HMA can yield averages of a given precision an order of magnitude
faster than conventional methods, and this only improves as the temperatures is
lowered.  This is particularly important for evaluating the free energy by
thermodynamic integration, where the low-temperature contributions are the
greatest source of statistical uncertainty.  Moreover, HMA has other
advantages, including smaller potential-truncation effects, finite-size
effects, smaller timestep inaccuracy, faster equilibration and shorter
decorrelation time.

HMA should not be used if atoms are expected to diffuse.  It is also
restricted to simulations in the NVT ensemble.  While this compute may be
used with any potential in LAMMPS, it will provide inaccurate results
for potentials that do not go to 0 at the truncation distance;
:doc:`pair\_lj\_smooth\_linear <pair_lj_smooth_linear>` and Ewald summation should
work fine, while :doc:`pair\_lj <pair_lj>` will perform poorly unless 
the potential is shifted (via :doc:`pair\_modify <pair_modify>` shift) or the cutoff is large.  Furthermore, computation of the heat capacity with
this compute is restricted to those that implement the single\_hessian method
in Pair.  Implementing single\_hessian in additional pair styles is simple.
Please contact Andrew Schultz (ajs42 at buffalo.edu) and David Kofke (kofke at
buffalo.edu) if your desired pair style does not have this method.  This is
the list of pair styles that currently implement pair\_hessian:

* :doc:`lj\_smooth\_linear <pair_lj_smooth_linear>`


In this method, the analytically known harmonic behavior of a crystal is removed from the traditional ensemble
averages, which leads to an accurate and precise measurement of the anharmonic contributions without contamination 
by noise produced by the already-known harmonic behavior. 
A detailed description of this method can be found in (:ref:`Moustafa <hma-Moustafa>`). The potential energy is computed by the formula:


.. math::

   \begin{equation}\left< U\right>_{HMA} = \frac{d}{2} (N-1) k_B T  + \left< U + \frac{1}{2} F\bullet\Delta r \right>\end{equation}

where :math:`N` is the number of atoms in the system, :math:`k_B` is Boltzmann's
constant, :math:`T` is the temperature, :math:`d` is the
dimensionality of the system (2 or 3 for 2d/3d), :math:`F\bullet\Delta r` is the sum of dot products of the 
atomic force vectors and displacement (from lattice sites) vectors, and :math:`U` is the sum of 
pair, bond, angle, dihedral, improper, kspace (long-range), and fix energies.

The pressure is computed by the formula:


.. math::

   \begin{equation}\left< P\right>_{HMA} = \Delta \hat P + \left< P_{vir} + \frac{\beta \Delta \hat P - \rho}{d(N-1)} F\bullet\Delta r \right>\end{equation}

where :math:`\rho` is the number density of the system, :math:`\Delta \hat P` is the
difference between the harmonic and lattice pressure, :math:`P_{vir}` is
the virial pressure computed as the sum of pair, bond, angle, dihedral,
improper, kspace (long-range), and fix contributions to the force on each
atom, and :math:`k_B=1/k_B T`.  Although the method will work for any value of :math:`\Delta \hat P`
specified (use pressure :doc:`units <units>`), the precision of the resultant
pressure is sensitive to :math:`\Delta \hat P`; the precision tends to be
best when :math:`\Delta \hat P` is the actual the difference between the lattice
pressure and harmonic pressure.


.. math::

   \begin{equation}\left<C_V \right>_{HMA} = \frac{d}{2} (N-1) k_B + \frac{1}{k_B T^2} \left( \left<
   U_{HMA}^2 \right> - \left<U_{HMA}\right>^2 \right) + \frac{1}{4 T}
   \left< F\bullet\Delta r + \Delta r \bullet \Phi \bullet \Delta r \right>\end{equation}

where :math:`\Phi` is the Hessian matrix. The compute hma command
computes the full expression for :math:`C_V` except for the
:math:`\left<U_{HMA}^2\right>^2` in the variance term, which can be obtained by
passing the *u* keyword; you must add this extra contribution to the :math:`C_V`
value reported by this compute.  The variance term can cause significant
round-off error when computing :math:`C_V`.  To address this, the *anharmonic*
keyword can be passed and/or the output format can be specified with more
digits.


.. parsed-literal::

   thermo_modify format float '%22.15e'

The *anharmonic* keyword will instruct the compute to return anharmonic
properties rather than the full properties, which include lattice, harmonic
and anharmonic contributions.
When using this keyword, the compute must be first active (it must be included
via a :doc:`thermo\_style custom <thermo_style>` command) while the atoms are
still at their lattice sites (before equilibration).

The temp-ID specified with compute hma command should be same as the fix-ID of Nose-Hoover (:doc:`fix nvt <fix_nh>`) or 
Berendsen (:doc:`fix temp/berendsen <fix_temp_berendsen>`) thermostat used for the simulation. While using this command, Langevin thermostat 
(:doc:`fix langevin <fix_langevin>`) 
should be avoided as its extra forces interfere with the HMA implementation.

.. note::

   Compute hma command should be used right after the energy minimization, when the atoms are at their lattice sites. 
   The simulation should not be started before this command has been used in the input script.

The following example illustrates the placement of this command in the input script:


.. parsed-literal::

   min_style cg 
   minimize 1e-35 1e-15 50000 500000 
   compute 1 all hma thermostatid u
   fix thermostatid all nvt temp 600.0 600.0 100.0

.. note::

   Compute hma should be used when the atoms of the solid do not diffuse. Diffusion will reduce the precision in the potential energy computation.

.. note::

   The :doc:`fix\_modify energy yes <fix_modify>` command must also be specified if a fix is to contribute potential energy to this command.

An example input script that uses this compute is included in
examples/USER/hma/ along with corresponding LAMMPS output showing that the HMA
properties fluctuate less than the corresponding conventional properties.

**Output info:**

This compute calculates a global vector that includes the n properties
requested as arguments to the command (the potential energy, pressure and/or heat
capacity).  The elements of the vector can be accessed by indices 1-n by any
command that uses global vector values as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output options.

The vector values calculated by this compute are "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-MISC package.  It is enabled only
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Usage restricted to canonical (NVT) ensemble simulation only.

Related commands
""""""""""""""""

:doc:`compute pe <compute_pe>`, :doc:`compute pressure <compute_pressure>`

:doc:`dynamical matrix <dynamical_matrix>` provides a finite difference
formulation of the hessian provided by Pair's single\_hessian, which is used by
this compute.

**Default:** none


----------


.. _hma-Moustafa:



**(Moustafa)** Sabry G. Moustafa, Andrew J. Schultz, and David A. Kofke, *Very fast averaging of thermal properties of crystals by molecular simulation*\ , 
`Phys. Rev. E [92], 043303 (2015) <https://link.aps.org/doi/10.1103/PhysRevE.92.043303>`_


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
