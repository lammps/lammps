.. table_from_list::
   :columns: 3

   * :doc:`General commands <Commands_all>`
   * :doc:`Fix styles <Commands_fix>`
   * :doc:`Compute styles <Commands_compute>`
   * :doc:`Pair styles <Commands_pair>`
   * :doc:`Bond styles <Commands_bond>`
   * :doc:`Angle styles <Commands_angle>`
   * :doc:`Dihedral styles <Commands_dihedral>`
   * :doc:`Improper styles <Commands_improper>`
   * :doc:`KSpace styles <Commands_kspace>`

KSpace solvers
==============

All LAMMPS :doc:`kspace\_style <kspace_style>` solvers.  Some styles have
accelerated versions.  This is indicated by additional letters in
parenthesis: g = GPU, i = USER-INTEL, k = KOKKOS, o = USER-OMP, t =
OPT.

.. table_from_list::
   :columns: 4

   * :doc:`ewald (o) <kspace_style>`
   * :doc:`ewald/disp <kspace_style>`
   * :doc:`msm (o) <kspace_style>`
   * :doc:`msm/cg (o) <kspace_style>`
   * :doc:`pppm (gok) <kspace_style>`
   * :doc:`pppm/cg (o) <kspace_style>`
   * :doc:`pppm/disp (i) <kspace_style>`
   * :doc:`pppm/disp/tip4p <kspace_style>`
   * :doc:`pppm/stagger <kspace_style>`
   * :doc:`pppm/tip4p (o) <kspace_style>`
   * :doc:`scafacos <kspace_style>`
   *
