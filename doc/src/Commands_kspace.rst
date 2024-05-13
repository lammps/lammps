.. table_from_list::
   :columns: 3

   * :doc:`General commands <Commands_all>`
   * :doc:`Fix styles <Commands_fix>`
   * :doc:`Compute styles <Commands_compute>`
   * :doc:`Pair styles <Commands_pair>`
   * :ref:`Bond styles <bond>`
   * :ref:`Angle styles <angle>`
   * :ref:`Dihedral styles <dihedral>`
   * :ref:`Improper styles <improper>`
   * :doc:`KSpace styles <Commands_kspace>`
   * :doc:`Dump styles <Commands_dump>`

KSpace solvers
==============

All LAMMPS :doc:`kspace_style <kspace_style>` solvers.  Some styles have
accelerated versions.  This is indicated by additional letters in
parenthesis: g = GPU, i = INTEL, k = KOKKOS, o = OPENMP, t =
OPT.

.. table_from_list::
   :columns: 4

   * :doc:`ewald (o) <kspace_style>`
   * :doc:`ewald/disp <kspace_style>`
   * :doc:`ewald/disp/dipole <kspace_style>`
   * :doc:`ewald/dipole <kspace_style>`
   * :doc:`ewald/dipole/spin <kspace_style>`
   * :doc:`ewald/electrode <kspace_style>`
   * :doc:`msm (o) <kspace_style>`
   * :doc:`msm/cg (o) <kspace_style>`
   * :doc:`msm/dielectric <kspace_style>`
   * :doc:`pppm (giko) <kspace_style>`
   * :doc:`pppm/cg (o) <kspace_style>`
   * :doc:`pppm/dipole <kspace_style>`
   * :doc:`pppm/dipole/spin <kspace_style>`
   * :doc:`pppm/dielectric <kspace_style>`
   * :doc:`pppm/disp (io) <kspace_style>`
   * :doc:`pppm/disp/tip4p (o) <kspace_style>`
   * :doc:`pppm/disp/dielectric <kspace_style>`
   * :doc:`pppm/stagger <kspace_style>`
   * :doc:`pppm/tip4p (o) <kspace_style>`
   * :doc:`pppm/dielectric <kspace_style>`
   * :doc:`pppm/electrode (i) <kspace_style>`
   * :doc:`scafacos <kspace_style>`
