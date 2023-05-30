.. index:: compute stress/cartesian
.. index:: compute stress/cylinder
.. index:: compute stress/spherical


compute stress/cartesian command
==================================

compute stress/cylinder command
=================================

compute stress/spherical command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID style args

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = stress/cartesian or stress/spherical or stress/cylinder
* args = argument specific to the compute style

.. parsed-literal::

  *stress/cartesian* args = dim bin_width
    dim = *x* or *y* or *z*
    bin_width = width of the bin
    one or two dim/bin_width pairs may be appended
  *stress/cylinder* args = zlo zh Rmax bin_width keyword
    zlo = minimum z-boundary for cylinder
    zhi = maximum z-boundary for cylinder
    Rmax = maximum radius to perform calculation to
    bin_width = width of radial bins to use for calculation
    keyword = ke (zero or one can be specified)
      ke = yes or no
  *stress/spherical*
    x0, y0, z0 = origin of the spherical coordinate system
    bin_width = width of spherical shells
    Rmax = maximum radius of spherical shells

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all stress/cartesian x 0.1
   compute 1 all stress/cartesian y 0.25 z 0.1
   compute 1 all stress/cylinder -10.0 10.0 15.0 0.25
   compute 1 all stress/cylinder -10.0 10.0 15.0 0.25 ke no
   compute 1 all stress/spherical 0 0 0 0.1 10

Description
"""""""""""

Compute *stress/cartesian*, compute *stress/cylinder*, and compute
*stress/spherical* define computations that calculate profiles of the
diagonal components of the local stress tensor in the specified
coordinate system. The stress tensor is split into a kinetic
contribution :math:`P^k` and a virial contribution :math:`P^v`. The sum
gives the total stress tensor :math:`P = P^k+P^v`. These computes can
for example be used to calculate the diagonal components of the local
stress tensor of interfaces with flat, cylindrical, or spherical
symmetry. These computes obeys momentum balance through fluid
interfaces. They use the Irving--Kirkwood contour, which is the straight
line between particle pairs.

The *stress/cartesian* computes the stress profile along one or two
Cartesian coordinates, as described in :ref:`(Ikeshoji)<Ikeshoji2>`. The
compute *stress/cylinder* computes the stress profile along the
radial direction in cylindrical coordinates, as described in
:ref:`(Addington)<Addington1>`. The compute *stress/spherical*
computes the stress profile along the radial direction in spherical
coordinates, as described in :ref:`(Ikeshoji)<Ikeshoji2>`.


Output info
"""""""""""

The output columns for *stress/cartesian* are the position of the
center of the local volume in the first and second dimensions, number
density, :math:`P^k_{xx}`, :math:`P^k_{yy}`, :math:`P^k_{zz}`,
:math:`P^v_{xx}`, :math:`P^v_{yy}`, and :math:`P^v_{zz}`. There are 8
columns when one dimension is specified and 9 columns when two
dimensions are specified. The number of bins (rows) is
:math:`(L_1/b_1)(L_2/b_2)`, where :math:`L_1` and :math:`L_2` are the lengths
of the simulation box in the specified dimensions and :math:`b_1` and
:math:`b_2` are the specified bin widths. When only one dimension is
specified, the number of bins (rows) is :math:`L_1/b_1`.

The default output columns for *stress/cylinder* are the radius to the
center of the cylindrical shell, number density, :math:`P^k_{rr}`,
:math:`P^k_{\phi\phi}`, :math:`P^k_{zz}`, :math:`P^v_{rr}`,
:math:`P^v_{\phi\phi}`, and :math:`P^v_{zz}`. When the keyword *ke* is
set to *no*, the kinetic contributions are not calculated, and
consequently there are only 5 columns: the position of the center of the
cylindrical shell, the number density, :math:`P^v_{rr}`,
:math:`P^v_{\phi\phi}`, and :math:`P^v_{zz}`. The number of bins (rows) is
:math:`R_\text{max}/b`, where :math:`b` is the specified bin width.

The output columns for *stress/spherical* are the position of the center
of the spherical shell, the number density, :math:`P^k_{rr}`,
:math:`P^k_{\theta\theta}`, :math:`P^k_{\phi\phi}`, :math:`P^v_{rr}`,
:math:`P^v_{\theta\theta}`, and :math:`P^v_{\phi\phi}`. There are 8
columns and the number of bins (rows) is :math:`R_\text{max}/b`, where
:math:`b` is the specified bin width.

This array can be output with :doc:`fix ave/time <fix_ave_time>`,

.. code-block:: LAMMPS

  compute p all stress/cartesian x 0.1
  fix 2 all ave/time 100 1 100 c_p[*] file dump_p.out mode vector

The values calculated by this compute are "intensive".  The stress
values will be in pressure :doc:`units <units>`. The number density
values are in inverse volume :doc:`units <units>`.

NOTE 1: The local stress does not include any Lennard-Jones tail
corrections to the stress added by the
:doc:`pair_modify tail yes <pair_modify>`
command, since those are contributions to the global system pressure.

NOTE 2: The local stress profiles generated by these computes are
similar to those obtained by the
:doc:`method-of-planes (MOP) <compute_stress_mop>`.
A key difference
is that compute `stress/mop/profile <compute_stress_mop>`
considers particles crossing a set of planes, while
*stress/cartesian* computes averages for a set of small volumes.
More information on the similarities and differences can be found in
:ref:`(Ikeshoji)<Ikeshoji2>`.

Restrictions
""""""""""""

These computes calculate the stress tensor contributions for pair styles
only (i.e., no bond, angle, dihedral, etc. contributions, and in the
presence of bonded interactions, the result may be incorrect due to
exclusions for :doc:`special bonds <special_bonds>` excluding pairs of atoms
completely). It requires pairwise force calculations not available for most
many-body pair styles.  Note that :math:`k`-space calculations are also excluded.

These computes are part of the EXTRA-COMPUTE package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build
package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`, :doc:`compute pressure <compute_pressure>`, :doc:`compute stress/mop/profile <compute_stress_mop>`

Default
"""""""

The keyword default for ke in style *stress/cylinder* is yes.

----------

.. _Ikeshoji2:

**(Ikeshoji)** Ikeshoji, Hafskjold, Furuholt, Mol Sim, 29, 101-109, (2003).

.. _Addington1:

**(Addington)** Addington, Long, Gubbins, J Chem Phys, 149, 084109 (2018).
