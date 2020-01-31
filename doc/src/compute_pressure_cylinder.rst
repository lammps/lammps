.. index:: compute pressure/cylinder

compute pressure/cylinder command
=================================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID pressure/cylinder zlo zhi Rmax bin_width

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/cylinder = style name of this compute command
* zlo = minimum z-boundary for cylinder
* zhi = maximum z-boundary for cylinder
* Rmax = maximum radius to perform calculation to
* bin\_width = width of radial bins to use for calculation

Examples
""""""""


.. parsed-literal::

   compute 1 all pressure/cylinder -10.0 10.0 15.0 0.25

Description
"""""""""""

Define a computation that calculates the pressure tensor of a system in
cylindrical coordinates, as discussed in :ref:`(Addington) <Addington1>`.
This is useful for systems with a single axis of rotational symmetry,
such as cylindrical micelles or carbon nanotubes. The compute splits the
system into radial, cylindrical-shell-type bins of width bin\_width,
centered at x=0,y=0, and calculates the radial (P\_rhorho), azimuthal
(P\_phiphi), and axial (P\_zz) components of the configurational pressure
tensor. The local density is also calculated for each bin, so that the
true pressure can be recovered as P\_kin+P\_conf=density\*k\*T+P\_conf.  The
output is a global array with 5 columns; one each for bin radius, local
number density, P\_rhorho, P\_phiphi, and P\_zz. The number of rows is
governed by the values of Rmax and bin\_width. Pressure tensor values are
output in pressure units.

**Output info:**

This compute calculates a global array with 5 columns and Rmax/bin\_width
rows. The output columns are: R (distance units), number density (inverse
volume units), configurational radial pressure (pressure units),
configurational azimuthal pressure (pressure units), and configurational
axial pressure (pressure units).

The values calculated by this compute are
"intensive".  The pressure values will be in pressure
:doc:`units <units>`. The number density values will be in
inverse volume :doc:`units <units>`.

Restrictions
""""""""""""


This compute currently calculates the pressure tensor contributions
for pair styles only (i.e. no bond, angle, dihedral, etc. contributions
and in the presence of bonded interactions, the result will be incorrect
due to exclusions for special bonds)  and requires pair-wise force
calculations not available for most many-body pair styles. K-space
calculations are also excluded. Note that this pressure compute outputs
the configurational terms only; the kinetic contribution is not included
and may be calculated from the number density output by P\_kin=density\*k\*T.

This compute is part of the USER-MISC package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute stress/atom <compute_stress_atom>`,
:doc:`thermo_style <thermo_style>`,

**Default:** none


----------


.. _Addington1:



**(Addington)** Addington, Long, Gubbins, J Chem Phys, 149, 084109 (2018).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
