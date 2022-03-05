.. index:: compute pressure/cartesian
.. index:: compute pressure/cylinder
.. index:: compute pressure/spherical


compute pressure/cartesian command
==================================

compute pressure/cylinder command
=================================

compute pressure/spherical command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID style args

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = pressure/cartesian or pressure/spherical or pressure/cylinder
* args = argument specific to the compute style

.. parsed-literal::

  *pressure/cartesian* args = dim bin_width
    dim = x, y, or z. One or two dim/bin_width pairs may be appended
    bin_width = width of the bin
  *pressure/cylinder* args = zlo zh Rmax bin_width keyword
    zlo = minimum z-boundary for cylinder
    zhi = maximum z-boundary for cylinder
    Rmax = maximum radius to perform calculation to
    bin_width = width of radial bins to use for calculation
    keyword = ke (zero or one can be specified)
      ke = yes or no
  *pressure/spherical*
    x0, y0, z0 = origin of the spherical coordinate system
    bin_width = width of spherical shells
    Rmax = maximum radius of spherical shells

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all pressure/cartesian x 0.1
   compute 1 all pressure/cartesian y 0.25 z 0.1
   compute 1 all pressure/cylinder -10.0 10.0 15.0 0.25
   compute 1 all pressure/cylinder -10.0 10.0 15.0 0.25 ke no
   compute 1 all pressure/spherical 0 0 0 0.1 10

Description
"""""""""""

Compute *pressure/cartesian*, compute *pressure/cylinder*, and compute *pressure/spherical* define computations that calculate profiles of the diagonal components of the local pressure tensor in the specified coordiante system. The pressure tensor is split into a kinetic contribution :math:`P^k` and a virial contribution :math:`P^v`. The sum gives the total pressure tensor :math:`P = P^k+P^v`. These computes can for example be used to calculate the diagonal components of the local pressure tensor of interfaces with flat, cylindrical, or spherical symmetry. These computes obeys momentum balance through fluid interfaces. These computes uses the Irving-Kirkwood contour, which is the straight line between particle pairs. 

The *pressure/cartesian* computes the pressure profile along one or two Cartesian corodinates, as described in :ref:`(Ikeshoji)<Ikeshoji2>`. The compute *pressure/cylinder* computes the pressure profile along the radial direction in cylindrical coordinates, as described in :ref:`(Addington)<Addington1>`. The compute *pressure/spherical* computes the pressure profile along the radial direction in spherical coordinates, as described in :ref:`(Ikeshoji)<Ikeshoji2>`.


Output info
"""""""""""

The output columns for *pressure/cartesian* are: position of the center of the local volume in the first and second dimensions, local number density, :math:`P^k_{xx}`, :math:`P^k_{yy}`, :math:`P^k_{zz}`, :math:`P^v_{xx}`, :math:`P^v_{yy}`, and :math:`P^v_{zz}`. There are 8 columns when one dimension is specified and 9 columns when two dimensions are specified. The number of bins/rows are (L1/bin_width1)*(L2/bin_width2), L1 and L2 are the size of the simulation box in the specified dimensions, and bin_widt1 and bin_width2 are the specified bin widths. When only one dimension is specified the number of bins/rows are L1/bin_width.

The default output columns for *pressure/cylinder* are: position of the center of the local volume, local number density, :math:`P^k_{rr}`, :math:`P^k_{\phi\phi}`, :math:`P^k_{zz}`, :math:`P^v_{rr}`, :math:`P^v_{\phi\phi}`, and :math:`P^v_{zz}`. When the keyword *ke* is set to no, the kinetic contributions are not calculated and consequently there are only 5 columns: radius to the center of the cylindrical shell, number density, :math:`P^v_{rr}`, :math:`P^v_{\phi\phi}`, :math:`P^v_{zz}`. The number of bins/rows are Rmax/bin_width.

The output columns for *pressure/spherical* are: radius to the center of the spherical shell, number density, :math:`P^k_{rr}`, :math:`P^k_{\theta\theta}`, :math:`P^k_{\phi\phi}`, :math:`P^v_{rr}`, :math:`P^v_{\theta\theta}`, and :math:`P^v_{\phi\phi}`. There are 8 columns and the number of bins/rows are Rmax/bin_width.

This array can be output with :doc:`fix ave/time <fix_ave_time>`,

.. code-block:: LAMMPS

  compute p all pressure/cartesian x 0.1
  fix 2 all ave/time 100 1 100 c_p[*] file dump_p.out mode vector

The values calculated by this compute are "intensive".  The pressure values will be in pressure :doc:`units <units>`. The number density values are in inverse volume :doc:`units <units>`.


Restrictions
""""""""""""

This compute currently calculates the pressure tensor contributions for pair styles only (i.e. no bond, angle, dihedral, etc. contributions and in the presence of bonded interactions, the result will be incorrect due to exclusions for special bonds) and requires pairwise force calculations not available for most many-body pair styles. K-space calculations are also excluded.

This compute is part of the EXTRA-COMPUTE package.  It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`, :doc:`compute pressure <compute_pressure>`, :doc:`compute stress/mop <compute_stress_mop>`

Default
"""""""

The keyword default for ke in style *pressure/cylinder* is yes.

----------

.. _Ikeshoji2:

**(Ikeshoji)** Ikeshoji, Hafskjold, Furuholt, Mol Sim, 29, 101-109, (2003).

.. _Addington1:

**(Addington)** Addington, Long, Gubbins, J Chem Phys, 149, 084109 (2018).
