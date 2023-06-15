.. index:: compute stress/cylinder
.. index:: compute stress/spherical

compute stress/cylinder command
=================================

compute stress/spherical command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID style args

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = stress/spherical or stress/cylinder
* args = argument specific to the compute style

.. parsed-literal::

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

   compute 1 all stress/cylinder -10.0 10.0 15.0 0.25
   compute 1 all stress/cylinder -10.0 10.0 15.0 0.25 ke no
   compute 1 all stress/spherical 0 0 0 0.1 10

Description
"""""""""""

Compute *stress/cylinder*, and compute
*stress/spherical* define computations that calculate profiles of the
diagonal components of the local stress tensor in the specified
coordinate system. The stress tensor is split into a kinetic
contribution :math:`P^k` and a virial contribution :math:`P^v`. The sum
gives the total stress tensor :math:`P = P^k+P^v`. These computes can
for example be used to calculate the diagonal components of the local
stress tensor of surfaces with cylindrical or spherical
symmetry. These computes obeys momentum balance through fluid
interfaces. They use the Irving--Kirkwood contour, which is the straight
line between particle pairs.

The compute *stress/cylinder* computes the stress profile along the
radial direction in cylindrical coordinates, as described in
:ref:`(Addington)<Addington1>`. The compute *stress/spherical*
computes the stress profile along the radial direction in spherical
coordinates, as described in :ref:`(Ikeshoji)<Ikeshoji4>`.


Output info
"""""""""""

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

  compute 1 all stress/spherical 0 0 0 0.1 10
  fix 2 all ave/time 100 1 100 c_p[*] file dump_p.out mode vector

The values calculated by this compute are "intensive".  The stress
values will be in pressure :doc:`units <units>`. The number density
values are in inverse volume :doc:`units <units>`.

NOTE 1: The local stress does not include any Lennard-Jones tail
corrections to the stress added by the
:doc:`pair_modify tail yes <pair_modify>`
command, since those are contributions to the global system pressure.

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

:doc:`compute stress/atom <compute_stress_atom>`, :doc:`compute pressure <compute_pressure>`,
:doc:`compute stress/mop/profile <compute_stress_mop>`, :doc:`compute stress/cartesian <compute_stress_cartesian>`

Default
"""""""

The keyword default for ke in style *stress/cylinder* is yes.

----------

.. _Ikeshoji4:

**(Ikeshoji)** Ikeshoji, Hafskjold, Furuholt, Mol Sim, 29, 101-109, (2003).

.. _Addington1:

**(Addington)** Addington, Long, Gubbins, J Chem Phys, 149, 084109 (2018).
