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

Fix commands
============

An alphabetic list of all LAMMPS :doc:`fix <fix>` commands.  Some styles
have accelerated versions.  This is indicated by additional letters in
parenthesis: g = GPU, i = INTEL, k = KOKKOS, o = OPENMP, t =
OPT.

.. table_from_list::
   :columns: 5

   * :doc:`accelerate/cos <fix_accelerate_cos>`
   * :doc:`acks2/reaxff (k) <fix_acks2_reaxff>`
   * :doc:`adapt <fix_adapt>`
   * :doc:`adapt/fep <fix_adapt_fep>`
   * :doc:`addforce <fix_addforce>`
   * :doc:`addtorque <fix_addtorque>`
   * :doc:`alchemy <fix_alchemy>`
   * :doc:`amoeba/bitorsion <fix_amoeba_bitorsion>`
   * :doc:`amoeba/pitorsion <fix_amoeba_pitorsion>`
   * :doc:`append/atoms <fix_append_atoms>`
   * :doc:`atc <fix_atc>`
   * :doc:`atom/swap <fix_atom_swap>`
   * :doc:`ave/atom <fix_ave_atom>`
   * :doc:`ave/chunk <fix_ave_chunk>`
   * :doc:`ave/correlate <fix_ave_correlate>`
   * :doc:`ave/correlate/long <fix_ave_correlate_long>`
   * :doc:`ave/grid <fix_ave_grid>`
   * :doc:`ave/histo <fix_ave_histo>`
   * :doc:`ave/histo/weight <fix_ave_histo>`
   * :doc:`ave/time <fix_ave_time>`
   * :doc:`aveforce <fix_aveforce>`
   * :doc:`balance <fix_balance>`
   * :doc:`bocs <fix_bocs>`
   * :doc:`bond/break <fix_bond_break>`
   * :doc:`bond/create <fix_bond_create>`
   * :doc:`bond/create/angle <fix_bond_create>`
   * :doc:`bond/react <fix_bond_react>`
   * :doc:`bond/swap <fix_bond_swap>`
   * :doc:`box/relax <fix_box_relax>`
   * :doc:`brownian <fix_brownian>`
   * :doc:`brownian/asphere <fix_brownian>`
   * :doc:`brownian/sphere <fix_brownian>`
   * :doc:`charge/regulation <fix_charge_regulation>`
   * :doc:`cmap <fix_cmap>`
   * :doc:`colvars <fix_colvars>`
   * :doc:`controller <fix_controller>`
   * :doc:`damping/cundall <fix_damping_cundall>`
   * :doc:`deform (k) <fix_deform>`
   * :doc:`deposit <fix_deposit>`
   * :doc:`dpd/energy (k) <fix_dpd_energy>`
   * :doc:`drag <fix_drag>`
   * :doc:`drude <fix_drude>`
   * :doc:`drude/transform/direct <fix_drude_transform>`
   * :doc:`drude/transform/inverse <fix_drude_transform>`
   * :doc:`dt/reset (k) <fix_dt_reset>`
   * :doc:`edpd/source <fix_dpd_source>`
   * :doc:`efield <fix_efield>`
   * :doc:`efield/tip4p <fix_efield>`
   * :doc:`ehex <fix_ehex>`
   * :doc:`electrode/conp (i) <fix_electrode>`
   * :doc:`electrode/conq (i) <fix_electrode>`
   * :doc:`electrode/thermo (i) <fix_electrode>`
   * :doc:`electron/stopping <fix_electron_stopping>`
   * :doc:`electron/stopping/fit <fix_electron_stopping>`
   * :doc:`enforce2d (k) <fix_enforce2d>`
   * :doc:`eos/cv <fix_eos_cv>`
   * :doc:`eos/table <fix_eos_table>`
   * :doc:`eos/table/rx (k) <fix_eos_table_rx>`
   * :doc:`evaporate <fix_evaporate>`
   * :doc:`external <fix_external>`
   * :doc:`ffl <fix_ffl>`
   * :doc:`filter/corotate <fix_filter_corotate>`
   * :doc:`flow/gauss <fix_flow_gauss>`
   * :doc:`freeze (k) <fix_freeze>`
   * :doc:`gcmc <fix_gcmc>`
   * :doc:`gld <fix_gld>`
   * :doc:`gle <fix_gle>`
   * :doc:`gravity (ko) <fix_gravity>`
   * :doc:`grem <fix_grem>`
   * :doc:`halt <fix_halt>`
   * :doc:`heat <fix_heat>`
   * :doc:`heat/flow <fix_heat_flow>`
   * :doc:`hyper/global <fix_hyper_global>`
   * :doc:`hyper/local <fix_hyper_local>`
   * :doc:`imd <fix_imd>`
   * :doc:`indent <fix_indent>`
   * :doc:`ipi <fix_ipi>`
   * :doc:`langevin (k) <fix_langevin>`
   * :doc:`langevin/drude <fix_langevin_drude>`
   * :doc:`langevin/eff <fix_langevin_eff>`
   * :doc:`langevin/spin <fix_langevin_spin>`
   * :doc:`lb/fluid <fix_lb_fluid>`
   * :doc:`lb/momentum <fix_lb_momentum>`
   * :doc:`lb/viscous <fix_lb_viscous>`
   * :doc:`lineforce <fix_lineforce>`
   * :doc:`manifoldforce <fix_manifoldforce>`
   * :doc:`mdi/qm <fix_mdi_qm>`
   * :doc:`mdi/qmmm <fix_mdi_qmmm>`
   * :doc:`meso/move <fix_meso_move>`
   * :doc:`mol/swap <fix_mol_swap>`
   * :doc:`momentum (k) <fix_momentum>`
   * :doc:`momentum/chunk <fix_momentum>`
   * :doc:`move <fix_move>`
   * :doc:`msst <fix_msst>`
   * :doc:`mvv/dpd <fix_mvv_dpd>`
   * :doc:`mvv/edpd <fix_mvv_dpd>`
   * :doc:`mvv/tdpd <fix_mvv_dpd>`
   * :doc:`neb <fix_neb>`
   * :doc:`neb/spin <fix_neb_spin>`
   * :doc:`nph (ko) <fix_nh>`
   * :doc:`nph/asphere (o) <fix_nph_asphere>`
   * :doc:`nph/body <fix_nph_body>`
   * :doc:`nph/eff <fix_nh_eff>`
   * :doc:`nph/sphere (o) <fix_nph_sphere>`
   * :doc:`nphug <fix_nphug>`
   * :doc:`npt (giko) <fix_nh>`
   * :doc:`npt/asphere (o) <fix_npt_asphere>`
   * :doc:`npt/body <fix_npt_body>`
   * :doc:`npt/cauchy <fix_npt_cauchy>`
   * :doc:`npt/eff <fix_nh_eff>`
   * :doc:`npt/sphere (o) <fix_npt_sphere>`
   * :doc:`npt/uef <fix_nh_uef>`
   * :doc:`numdiff <fix_numdiff>`
   * :doc:`numdiff/virial <fix_numdiff_virial>`
   * :doc:`nve (giko) <fix_nve>`
   * :doc:`nve/asphere (gi) <fix_nve_asphere>`
   * :doc:`nve/asphere/noforce <fix_nve_asphere_noforce>`
   * :doc:`nve/awpmd <fix_nve_awpmd>`
   * :doc:`nve/body <fix_nve_body>`
   * :doc:`nve/dot <fix_nve_dot>`
   * :doc:`nve/dotc/langevin <fix_nve_dotc_langevin>`
   * :doc:`nve/eff <fix_nve_eff>`
   * :doc:`nve/limit <fix_nve_limit>`
   * :doc:`nve/line <fix_nve_line>`
   * :doc:`nve/manifold/rattle <fix_nve_manifold_rattle>`
   * :doc:`nve/noforce <fix_nve_noforce>`
   * :doc:`nve/sphere (ko) <fix_nve_sphere>`
   * :doc:`nve/bpm/sphere <fix_nve_bpm_sphere>`
   * :doc:`nve/spin <fix_nve_spin>`
   * :doc:`nve/tri <fix_nve_tri>`
   * :doc:`nvk <fix_nvk>`
   * :doc:`nvt (giko) <fix_nh>`
   * :doc:`nvt/asphere (o) <fix_nvt_asphere>`
   * :doc:`nvt/body <fix_nvt_body>`
   * :doc:`nvt/eff <fix_nh_eff>`
   * :doc:`nvt/manifold/rattle <fix_nvt_manifold_rattle>`
   * :doc:`nvt/sllod (iko) <fix_nvt_sllod>`
   * :doc:`nvt/sllod/eff <fix_nvt_sllod_eff>`
   * :doc:`nvt/sphere (o) <fix_nvt_sphere>`
   * :doc:`nvt/uef <fix_nh_uef>`
   * :doc:`oneway <fix_oneway>`
   * :doc:`orient/bcc <fix_orient>`
   * :doc:`orient/fcc <fix_orient>`
   * :doc:`orient/eco <fix_orient_eco>`
   * :doc:`pafi <fix_pafi>`
   * :doc:`pair <fix_pair>`
   * :doc:`phonon <fix_phonon>`
   * :doc:`pimd/langevin <fix_pimd>`
   * :doc:`pimd/nvt <fix_pimd>`
   * :doc:`planeforce <fix_planeforce>`
   * :doc:`plumed <fix_plumed>`
   * :doc:`poems <fix_poems>`
   * :doc:`polarize/bem/gmres <fix_polarize>`
   * :doc:`polarize/bem/icc <fix_polarize>`
   * :doc:`polarize/functional <fix_polarize>`
   * :doc:`pour <fix_pour>`
   * :doc:`precession/spin <fix_precession_spin>`
   * :doc:`press/berendsen <fix_press_berendsen>`
   * :doc:`print <fix_print>`
   * :doc:`propel/self <fix_propel_self>`
   * :doc:`property/atom (k) <fix_property_atom>`
   * :doc:`python/invoke <fix_python_invoke>`
   * :doc:`python/move <fix_python_move>`
   * :doc:`qbmsst <fix_qbmsst>`
   * :doc:`qeq/comb (o) <fix_qeq_comb>`
   * :doc:`qeq/dynamic <fix_qeq>`
   * :doc:`qeq/fire <fix_qeq>`
   * :doc:`qeq/point <fix_qeq>`
   * :doc:`qeq/reaxff (ko) <fix_qeq_reaxff>`
   * :doc:`qeq/shielded <fix_qeq>`
   * :doc:`qeq/slater <fix_qeq>`
   * :doc:`qmmm <fix_qmmm>`
   * :doc:`qtb <fix_qtb>`
   * :doc:`rattle <fix_shake>`
   * :doc:`reaxff/bonds (k) <fix_reaxff_bonds>`
   * :doc:`reaxff/species (k) <fix_reaxff_species>`
   * :doc:`recenter <fix_recenter>`
   * :doc:`restrain <fix_restrain>`
   * :doc:`rhok <fix_rhok>`
   * :doc:`rigid (o) <fix_rigid>`
   * :doc:`rigid/meso <fix_rigid_meso>`
   * :doc:`rigid/nph (o) <fix_rigid>`
   * :doc:`rigid/nph/small <fix_rigid>`
   * :doc:`rigid/npt (o) <fix_rigid>`
   * :doc:`rigid/npt/small <fix_rigid>`
   * :doc:`rigid/nve (o) <fix_rigid>`
   * :doc:`rigid/nve/small <fix_rigid>`
   * :doc:`rigid/nvt (o) <fix_rigid>`
   * :doc:`rigid/nvt/small <fix_rigid>`
   * :doc:`rigid/small (o) <fix_rigid>`
   * :doc:`rx (k) <fix_rx>`
   * :doc:`saed/vtk <fix_saed_vtk>`
   * :doc:`setforce (k) <fix_setforce>`
   * :doc:`setforce/spin <fix_setforce>`
   * :doc:`sgcmc <fix_sgcmc>`
   * :doc:`shake (k) <fix_shake>`
   * :doc:`shardlow (k) <fix_shardlow>`
   * :doc:`smd <fix_smd>`
   * :doc:`smd/adjust_dt <fix_smd_adjust_dt>`
   * :doc:`smd/integrate_tlsph <fix_smd_integrate_tlsph>`
   * :doc:`smd/integrate_ulsph <fix_smd_integrate_ulsph>`
   * :doc:`smd/move_tri_surf <fix_smd_move_triangulated_surface>`
   * :doc:`smd/setvel <fix_smd_setvel>`
   * :doc:`smd/wall_surface <fix_smd_wall_surface>`
   * :doc:`sph <fix_sph>`
   * :doc:`sph/stationary <fix_sph_stationary>`
   * :doc:`spring <fix_spring>`
   * :doc:`spring/chunk <fix_spring_chunk>`
   * :doc:`spring/rg <fix_spring_rg>`
   * :doc:`spring/self <fix_spring_self>`
   * :doc:`srd <fix_srd>`
   * :doc:`store/force <fix_store_force>`
   * :doc:`store/state <fix_store_state>`
   * :doc:`tdpd/source <fix_dpd_source>`
   * :doc:`temp/berendsen <fix_temp_berendsen>`
   * :doc:`temp/csld <fix_temp_csvr>`
   * :doc:`temp/csvr <fix_temp_csvr>`
   * :doc:`temp/rescale <fix_temp_rescale>`
   * :doc:`temp/rescale/eff <fix_temp_rescale_eff>`
   * :doc:`tfmc <fix_tfmc>`
   * :doc:`tgnpt/drude <fix_tgnh_drude>`
   * :doc:`tgnvt/drude <fix_tgnh_drude>`
   * :doc:`thermal/conductivity <fix_thermal_conductivity>`
   * :doc:`ti/spring <fix_ti_spring>`
   * :doc:`tmd <fix_tmd>`
   * :doc:`ttm <fix_ttm>`
   * :doc:`ttm/grid <fix_ttm>`
   * :doc:`ttm/mod <fix_ttm>`
   * :doc:`tune/kspace <fix_tune_kspace>`
   * :doc:`vector <fix_vector>`
   * :doc:`viscosity <fix_viscosity>`
   * :doc:`viscous (k) <fix_viscous>`
   * :doc:`viscous/sphere <fix_viscous_sphere>`
   * :doc:`wall/body/polygon <fix_wall_body_polygon>`
   * :doc:`wall/body/polyhedron <fix_wall_body_polyhedron>`
   * :doc:`wall/colloid <fix_wall>`
   * :doc:`wall/ees <fix_wall_ees>`
   * :doc:`wall/gran (k) <fix_wall_gran>`
   * :doc:`wall/gran/region <fix_wall_gran_region>`
   * :doc:`wall/harmonic <fix_wall>`
   * :doc:`wall/lj1043 <fix_wall>`
   * :doc:`wall/lj126 <fix_wall>`
   * :doc:`wall/lj93 (k) <fix_wall>`
   * :doc:`wall/lepton <fix_wall>`
   * :doc:`wall/morse <fix_wall>`
   * :doc:`wall/piston <fix_wall_piston>`
   * :doc:`wall/reflect (k) <fix_wall_reflect>`
   * :doc:`wall/reflect/stochastic <fix_wall_reflect_stochastic>`
   * :doc:`wall/region <fix_wall_region>`
   * :doc:`wall/region/ees <fix_wall_ees>`
   * :doc:`wall/srd <fix_wall_srd>`
   * :doc:`wall/table <fix_wall>`
   * :doc:`widom <fix_widom>`
