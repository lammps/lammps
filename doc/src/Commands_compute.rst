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

Compute commands
================

An alphabetic list of all LAMMPS :doc:`compute <compute>` commands.
Some styles have accelerated versions.  This is indicated by
additional letters in parenthesis: g = GPU, i = USER-INTEL, k =
KOKKOS, o = USER-OMP, t = OPT.

.. table_from_list::
   :columns: 5

   * :doc:`ackland/atom <compute_ackland_atom>`
   * :doc:`adf <compute_adf>`
   * :doc:`aggregate/atom <compute_cluster_atom>`
   * :doc:`angle <compute_angle>`
   * :doc:`angle/local <compute_angle_local>`
   * :doc:`angmom/chunk <compute_angmom_chunk>`
   * :doc:`basal/atom <compute_basal_atom>`
   * :doc:`body/local <compute_body_local>`
   * :doc:`bond <compute_bond>`
   * :doc:`bond/local <compute_bond_local>`
   * :doc:`centro/atom <compute_centro_atom>`
   * :doc:`centroid/stress/atom <compute_stress_atom>`
   * :doc:`chunk/atom <compute_chunk_atom>`
   * :doc:`chunk/spread/atom <compute_chunk_spread_atom>`
   * :doc:`cluster/atom <compute_cluster_atom>`
   * :doc:`cna/atom <compute_cna_atom>`
   * :doc:`cnp/atom <compute_cnp_atom>`
   * :doc:`com <compute_com>`
   * :doc:`com/chunk <compute_com_chunk>`
   * :doc:`contact/atom <compute_contact_atom>`
   * :doc:`coord/atom <compute_coord_atom>`
   * :doc:`damage/atom <compute_damage_atom>`
   * :doc:`dihedral <compute_dihedral>`
   * :doc:`dihedral/local <compute_dihedral_local>`
   * :doc:`dilatation/atom <compute_dilatation_atom>`
   * :doc:`dipole/chunk <compute_dipole_chunk>`
   * :doc:`displace/atom <compute_displace_atom>`
   * :doc:`dpd <compute_dpd>`
   * :doc:`dpd/atom <compute_dpd_atom>`
   * :doc:`edpd/temp/atom <compute_edpd_temp_atom>`
   * :doc:`entropy/atom <compute_entropy_atom>`
   * :doc:`erotate/asphere <compute_erotate_asphere>`
   * :doc:`erotate/rigid <compute_erotate_rigid>`
   * :doc:`erotate/sphere <compute_erotate_sphere>`
   * :doc:`erotate/sphere/atom <compute_erotate_sphere_atom>`
   * :doc:`event/displace <compute_event_displace>`
   * :doc:`fep <compute_fep>`
   * :doc:`force/tally <compute_tally>`
   * :doc:`fragment/atom <compute_cluster_atom>`
   * :doc:`global/atom <compute_global_atom>`
   * :doc:`group/group <compute_group_group>`
   * :doc:`gyration <compute_gyration>`
   * :doc:`gyration/chunk <compute_gyration_chunk>`
   * :doc:`gyration/shape <compute_gyration_shape>`
   * :doc:`gyration/shape/chunk <compute_gyration_shape_chunk>`
   * :doc:`heat/flux <compute_heat_flux>`
   * :doc:`heat/flux/tally <compute_tally>`
   * :doc:`hexorder/atom <compute_hexorder_atom>`
   * :doc:`hma <compute_hma>`
   * :doc:`improper <compute_improper>`
   * :doc:`improper/local <compute_improper_local>`
   * :doc:`inertia/chunk <compute_inertia_chunk>`
   * :doc:`ke <compute_ke>`
   * :doc:`ke/atom <compute_ke_atom>`
   * :doc:`ke/atom/eff <compute_ke_atom_eff>`
   * :doc:`ke/eff <compute_ke_eff>`
   * :doc:`ke/rigid <compute_ke_rigid>`
   * :doc:`meso/e/atom <compute_meso_e_atom>`
   * :doc:`meso/rho/atom <compute_meso_rho_atom>`
   * :doc:`meso/t/atom <compute_meso_t_atom>`
   * :doc:`momentum <compute_momentum>`
   * :doc:`msd <compute_msd>`
   * :doc:`msd/chunk <compute_msd_chunk>`
   * :doc:`msd/nongauss <compute_msd_nongauss>`
   * :doc:`omega/chunk <compute_omega_chunk>`
   * :doc:`orientorder/atom (k) <compute_orientorder_atom>`
   * :doc:`pair <compute_pair>`
   * :doc:`pair/local <compute_pair_local>`
   * :doc:`pe <compute_pe>`
   * :doc:`pe/atom <compute_pe_atom>`
   * :doc:`pe/mol/tally <compute_tally>`
   * :doc:`pe/tally <compute_tally>`
   * :doc:`plasticity/atom <compute_plasticity_atom>`
   * :doc:`pressure <compute_pressure>`
   * :doc:`pressure/cylinder <compute_pressure_cylinder>`
   * :doc:`pressure/uef <compute_pressure_uef>`
   * :doc:`property/atom <compute_property_atom>`
   * :doc:`property/chunk <compute_property_chunk>`
   * :doc:`property/local <compute_property_local>`
   * :doc:`ptm/atom <compute_ptm_atom>`
   * :doc:`rdf <compute_rdf>`
   * :doc:`reduce <compute_reduce>`
   * :doc:`reduce/chunk <compute_reduce_chunk>`
   * :doc:`reduce/region <compute_reduce>`
   * :doc:`rigid/local <compute_rigid_local>`
   * :doc:`saed <compute_saed>`
   * :doc:`slice <compute_slice>`
   * :doc:`smd/contact/radius <compute_smd_contact_radius>`
   * :doc:`smd/damage <compute_smd_damage>`
   * :doc:`smd/hourglass/error <compute_smd_hourglass_error>`
   * :doc:`smd/internal/energy <compute_smd_internal_energy>`
   * :doc:`smd/plastic/strain <compute_smd_plastic_strain>`
   * :doc:`smd/plastic/strain/rate <compute_smd_plastic_strain_rate>`
   * :doc:`smd/rho <compute_smd_rho>`
   * :doc:`smd/tlsph/defgrad <compute_smd_tlsph_defgrad>`
   * :doc:`smd/tlsph/dt <compute_smd_tlsph_dt>`
   * :doc:`smd/tlsph/num/neighs <compute_smd_tlsph_num_neighs>`
   * :doc:`smd/tlsph/shape <compute_smd_tlsph_shape>`
   * :doc:`smd/tlsph/strain <compute_smd_tlsph_strain>`
   * :doc:`smd/tlsph/strain/rate <compute_smd_tlsph_strain_rate>`
   * :doc:`smd/tlsph/stress <compute_smd_tlsph_stress>`
   * :doc:`smd/triangle/vertices <compute_smd_triangle_vertices>`
   * :doc:`smd/ulsph/num/neighs <compute_smd_ulsph_num_neighs>`
   * :doc:`smd/ulsph/strain <compute_smd_ulsph_strain>`
   * :doc:`smd/ulsph/strain/rate <compute_smd_ulsph_strain_rate>`
   * :doc:`smd/ulsph/stress <compute_smd_ulsph_stress>`
   * :doc:`smd/vol <compute_smd_vol>`
   * :doc:`snap <compute_sna_atom>`
   * :doc:`sna/atom <compute_sna_atom>`
   * :doc:`snad/atom <compute_sna_atom>`
   * :doc:`snav/atom <compute_sna_atom>`
   * :doc:`spin <compute_spin>`
   * :doc:`stress/atom <compute_stress_atom>`
   * :doc:`stress/mop <compute_stress_mop>`
   * :doc:`stress/mop/profile <compute_stress_mop>`
   * :doc:`stress/tally <compute_tally>`
   * :doc:`tdpd/cc/atom <compute_tdpd_cc_atom>`
   * :doc:`temp (k) <compute_temp>`
   * :doc:`temp/asphere <compute_temp_asphere>`
   * :doc:`temp/body <compute_temp_body>`
   * :doc:`temp/chunk <compute_temp_chunk>`
   * :doc:`temp/com <compute_temp_com>`
   * :doc:`temp/cs <compute_temp_cs>`
   * :doc:`temp/deform <compute_temp_deform>`
   * :doc:`temp/deform/eff <compute_temp_deform_eff>`
   * :doc:`temp/drude <compute_temp_drude>`
   * :doc:`temp/eff <compute_temp_eff>`
   * :doc:`temp/partial <compute_temp_partial>`
   * :doc:`temp/profile <compute_temp_profile>`
   * :doc:`temp/ramp <compute_temp_ramp>`
   * :doc:`temp/region <compute_temp_region>`
   * :doc:`temp/region/eff <compute_temp_region_eff>`
   * :doc:`temp/rotate <compute_temp_rotate>`
   * :doc:`temp/sphere <compute_temp_sphere>`
   * :doc:`temp/uef <compute_temp_uef>`
   * :doc:`ti <compute_ti>`
   * :doc:`torque/chunk <compute_torque_chunk>`
   * :doc:`vacf <compute_vacf>`
   * :doc:`vcm/chunk <compute_vcm_chunk>`
   * :doc:`viscosity/cos <compute_viscosity_cos>`
   * :doc:`voronoi/atom <compute_voronoi_atom>`
   * :doc:`xrd <compute_xrd>`
