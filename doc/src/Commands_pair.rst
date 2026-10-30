.. include:: Commands_html.rst

Pair styles
===========

All LAMMPS :doc:`pair_style <pair_style>` commands.  Some styles have
accelerated versions.  This is indicated by additional letters in
parenthesis: g = GPU, i = INTEL, k = KOKKOS, o = OPENMP, t = OPT.

.. table_from_list::
   :columns: 3

   * :doc:`none <pair_none>`
   * :doc:`zero <pair_zero>`
   * :doc:`zero/coul <pair_zero>`
   * :doc:`hybrid (ko) <pair_hybrid>`
   * :doc:`hybrid/molecular (o) <pair_hybrid>`
   * :doc:`hybrid/overlay (ko) <pair_hybrid>`
   * :doc:`hybrid/scaled (ko) <pair_hybrid>`
   * :doc:`kim <pair_kim>`
   * :doc:`lepton (o) <pair_lepton>`
   * :doc:`lepton/coul (o) <pair_lepton>`
   * :doc:`lepton/sphere (o) <pair_lepton>`
   * :doc:`list <pair_list>`
   * :doc:`python <pair_python>`
   * :doc:`table (gko) <pair_table>`
   * :doc:`tracker <pair_tracker>`
   *
   *
   *
   * :doc:`adp (ko) <pair_adp>`
   * :doc:`aip/water/2dm (t) <pair_aip_water_2dm>`
   * :doc:`airebo (io) <pair_airebo>`
   * :doc:`airebo/bc (o) <pair_airebo>`
   * :doc:`airebo/morse (io) <pair_airebo>`
   * :doc:`amoeba (g) <pair_amoeba>`
   * :doc:`atm (ko) <pair_atm>`
   * :doc:`beck (gko) <pair_beck>`
   * :doc:`body/nparticle (o) <pair_body_nparticle>`
   * :doc:`body/rounded/polygon (o) <pair_body_rounded_polygon>`
   * :doc:`body/rounded/polyhedron (o) <pair_body_rounded_polyhedron>`
   * :doc:`bondval (ko) <pair_bondval>`
   * :doc:`bondval/vec (ko) <pair_bondval>`
   * :doc:`bop <pair_bop>`
   * :doc:`born (gko) <pair_born>`
   * :doc:`born/coul/dsf (gko) <pair_born>`
   * :doc:`born/coul/dsf/cs (gko) <pair_cs>`
   * :doc:`born/coul/long (gko) <pair_born>`
   * :doc:`born/coul/long/cs (gko) <pair_cs>`
   * :doc:`born/coul/msm (o) <pair_born>`
   * :doc:`born/coul/wolf (gko) <pair_born>`
   * :doc:`born/coul/wolf/cs (gko) <pair_cs>`
   * :doc:`born/gauss (ko) <pair_born_gauss>`
   * :doc:`bpm/peri (o) <pair_bpm_peri>`
   * :doc:`bpm/spring (o) <pair_bpm_spring>`
   * :doc:`brownian (ko) <pair_brownian>`
   * :doc:`brownian/poly (o) <pair_brownian>`
   * :doc:`buck (giko) <pair_buck>`
   * :doc:`buck/coul/cut (giko) <pair_buck>`
   * :doc:`buck/coul/long (giko) <pair_buck>`
   * :doc:`buck/coul/long/cs (gko) <pair_cs>`
   * :doc:`buck/coul/msm (o) <pair_buck>`
   * :doc:`buck/long/coul/long (o) <pair_buck_long>`
   * :doc:`buck/mdf (ko) <pair_mdf>`
   * :doc:`buck6d/coul/gauss/dsf (ko) <pair_buck6d_coul_gauss>`
   * :doc:`buck6d/coul/gauss/long (ko) <pair_buck6d_coul_gauss>`
   * :doc:`colloid (gko) <pair_colloid>`
   * :doc:`comb (o) <pair_comb>`
   * :doc:`comb3 <pair_comb>`
   * :doc:`cosine/squared (ko) <pair_cosine_squared>`
   * :doc:`coul/ctip (ko) <pair_coul>`
   * :doc:`coul/cut (gko) <pair_coul>`
   * :doc:`coul/cut/dielectric (o) <pair_dielectric>`
   * :doc:`coul/cut/global (gko) <pair_coul>`
   * :doc:`coul/cut/soft (gko) <pair_fep_soft>`
   * :doc:`coul/cut/soft/gapsys (gko) <pair_fep_soft>`
   * :doc:`coul/debye (gko) <pair_coul>`
   * :doc:`coul/diel (gko) <pair_coul_diel>`
   * :doc:`coul/dsf (gko) <pair_coul>`
   * :doc:`coul/esp (ko) <pair_coul>`
   * :doc:`coul/exclude (o) <pair_coul>`
   * :doc:`coul/long (gko) <pair_coul>`
   * :doc:`coul/long/cs (gko) <pair_cs>`
   * :doc:`coul/long/dielectric (o) <pair_dielectric>`
   * :doc:`coul/long/soft (gko) <pair_fep_soft>`
   * :doc:`coul/msm (o) <pair_coul>`
   * :doc:`coul/slater/cut (gko) <pair_coul_slater>`
   * :doc:`coul/slater/long (gko) <pair_coul_slater>`
   * :doc:`coul/shield (ko) <pair_coul_shield>`
   * :doc:`coul/streitz (o) <pair_coul>`
   * :doc:`coul/tt (o) <pair_coul_tt>`
   * :doc:`coul/wolf (gko) <pair_coul>`
   * :doc:`coul/wolf/cs (gko) <pair_cs>`
   * :doc:`dispersion/d3 <pair_dispersion_d3>`
   * :doc:`dpd (giko) <pair_dpd>`
   * :doc:`dpd/coul/slater/long (gko) <pair_dpd_coul_slater_long>`
   * :doc:`dpd/ext (ko) <pair_dpd_ext>`
   * :doc:`dpd/ext/tstat (ko) <pair_dpd_ext>`
   * :doc:`dpd/fdt (o) <pair_dpd_fdt>`
   * :doc:`dpd/fdt/energy (k) <pair_dpd_fdt>`
   * :doc:`dpd/tstat (gko) <pair_dpd>`
   * :doc:`dsmc <pair_dsmc>`
   * :doc:`e3b (o) <pair_e3b>`
   * :doc:`drip (o) <pair_drip>`
   * :doc:`eam (gikot) <pair_eam>`
   * :doc:`eam/alloy (gikot) <pair_eam>`
   * :doc:`eam/cd <pair_eam>`
   * :doc:`eam/cd/old <pair_eam>`
   * :doc:`eam/fs (gikot) <pair_eam>`
   * :doc:`eam/fs/apip <pair_eam_apip>`
   * :doc:`eam/he (gkot) <pair_eam>`
   * :doc:`eam/apip <pair_eam_apip>`
   * :doc:`edip (ko) <pair_edip>`
   * :doc:`edip/multi (ko) <pair_edip>`
   * :doc:`edpd (go) <pair_mesodpd>`
   * :doc:`eff/cut <pair_eff>`
   * :doc:`eim (o) <pair_eim>`
   * :doc:`exp6/rx (k) <pair_exp6_rx>`
   * :doc:`extep (o) <pair_extep>`
   * :doc:`gauss (gko) <pair_gauss>`
   * :doc:`gauss/cut (ko) <pair_gauss>`
   * :doc:`gayberne (gio) <pair_gayberne>`
   * :doc:`gran/hertz/history (ko) <pair_gran>`
   * :doc:`gran/hooke (ko) <pair_gran>`
   * :doc:`gran/hooke/history (ko) <pair_gran>`
   * :doc:`granular <pair_granular>`
   * :doc:`granular/superellipsoid <pair_granular_superellipsoid>`
   * :doc:`gw (ko) <pair_gw>`
   * :doc:`gw/zbl (ko) <pair_gw>`
   * :doc:`harmonic/cut (ko) <pair_harmonic_cut>`
   * :doc:`hbond/dreiding/lj (o) <pair_hbond_dreiding>`
   * :doc:`hbond/dreiding/lj/angleoffset (o) <pair_hbond_dreiding>`
   * :doc:`hbond/dreiding/morse (o) <pair_hbond_dreiding>`
   * :doc:`hbond/dreiding/morse/angleoffset (o) <pair_hbond_dreiding>`
   * :doc:`hdnnp <pair_hdnnp>`
   * :doc:`hippo (g) <pair_amoeba>`
   * :doc:`ilp/graphene/hbn (t) <pair_ilp_graphene_hbn>`
   * :doc:`ilp/tmd (t) <pair_ilp_tmd>`
   * :doc:`kolmogorov/crespi/full <pair_kolmogorov_crespi_full>`
   * :doc:`kolmogorov/crespi/z (o) <pair_kolmogorov_crespi_z>`
   * :doc:`lambda/input/apip <pair_lambda_input_apip>`
   * :doc:`lambda/input/csp/apip <pair_lambda_input_apip>`
   * :doc:`lambda/zone/apip <pair_lambda_zone_apip>`
   * :doc:`lcbop (o) <pair_lcbop>`
   * :doc:`ldd <pair_ldd>`
   * :doc:`lebedeva/z (o) <pair_lebedeva_z>`
   * :doc:`lennard/mdf (ko) <pair_mdf>`
   * :doc:`line/lj (o) <pair_line_lj>`
   * :doc:`lj/charmm/coul/charmm (giko) <pair_charmm>`
   * :doc:`lj/charmm/coul/charmm/implicit (gko) <pair_charmm>`
   * :doc:`lj/charmm/coul/long (gikot) <pair_charmm>`
   * :doc:`lj/charmm/coul/long/soft (gko) <pair_fep_soft>`
   * :doc:`lj/charmm/coul/msm (o) <pair_charmm>`
   * :doc:`lj/charmmfsw/coul/charmmfsh (gko) <pair_charmm>`
   * :doc:`lj/charmmfsw/coul/long (gko) <pair_charmm>`
   * :doc:`lj/class2 (gko) <pair_class2>`
   * :doc:`lj/class2/coul/cut (gko) <pair_class2>`
   * :doc:`lj/class2/coul/cut/soft (ko) <pair_fep_soft>`
   * :doc:`lj/class2/coul/cut/soft/gapsys (gko) <pair_fep_soft>`
   * :doc:`lj/class2/coul/long (gko) <pair_class2>`
   * :doc:`lj/class2/coul/long/cs (gko) <pair_cs>`
   * :doc:`lj/class2/coul/long/soft (gko) <pair_fep_soft>`
   * :doc:`lj/class2/soft (gko) <pair_fep_soft>`
   * :doc:`lj/class2/soft/gapsys (gko) <pair_fep_soft>`
   * :doc:`lj/cubic (gko) <pair_lj_cubic>`
   * :doc:`lj/cut (gikot) <pair_lj>`
   * :doc:`lj/cut/coul/cut (gko) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/cut/dielectric (o) <pair_dielectric>`
   * :doc:`lj/cut/coul/cut/soft (gko) <pair_fep_soft>`
   * :doc:`lj/cut/coul/cut/soft/gapsys (gko) <pair_fep_soft>`
   * :doc:`lj/cut/coul/debye (gko) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/debye/dielectric (o) <pair_dielectric>`
   * :doc:`lj/cut/coul/dsf (gko) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/esp (ko) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/long (gikot) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/long/cs (gko) <pair_cs>`
   * :doc:`lj/cut/coul/long/dielectric (o) <pair_dielectric>`
   * :doc:`lj/cut/coul/long/soft (gko) <pair_fep_soft>`
   * :doc:`lj/cut/coul/msm (go) <pair_lj_cut_coul>`
   * :doc:`lj/cut/coul/msm/dielectric (o) <pair_dielectric>`
   * :doc:`lj/cut/coul/wolf (gko) <pair_lj_cut_coul>`
   * :doc:`lj/cut/dipole/cut (gko) <pair_dipole>`
   * :doc:`lj/cut/dipole/long (go) <pair_dipole>`
   * :doc:`lj/cut/dipole/sf (go) <pair_dipole>`
   * :doc:`lj/cut/soft (gko) <pair_fep_soft>`
   * :doc:`lj/cut/soft/gapsys (gko) <pair_fep_soft>`
   * :doc:`lj/cut/sphere (ko) <pair_lj_cut_sphere>`
   * :doc:`lj/cut/thole/long (o) <pair_thole>`
   * :doc:`lj/cut/tip4p/cut (ko) <pair_lj_cut_tip4p>`
   * :doc:`lj/cut/tip4p/long (gkot) <pair_lj_cut_tip4p>`
   * :doc:`lj/cut/tip4p/long/soft (ko) <pair_fep_soft>`
   * :doc:`lj/expand (gko) <pair_lj_expand>`
   * :doc:`lj/expand/coul/long (gko) <pair_lj_expand>`
   * :doc:`lj/expand/sphere (ko) <pair_lj_expand_sphere>`
   * :doc:`lj/gromacs (gko) <pair_gromacs>`
   * :doc:`lj/gromacs/coul/gromacs (gko) <pair_gromacs>`
   * :doc:`lj/long/coul/long (iot) <pair_lj_long>`
   * :doc:`lj/long/coul/long/dielectric (o) <pair_dielectric>`
   * :doc:`lj/long/dipole/long (o) <pair_dipole>`
   * :doc:`lj/long/tip4p/long (o) <pair_lj_long>`
   * :doc:`lj/mdf (ko) <pair_mdf>`
   * :doc:`lj/pirani (ko) <pair_lj_pirani>`
   * :doc:`lj/relres (ko) <pair_lj_relres>`
   * :doc:`lj/spica (gko) <pair_spica>`
   * :doc:`lj/spica/coul/long (gko) <pair_spica>`
   * :doc:`lj/spica/coul/msm (o) <pair_spica>`
   * :doc:`lj/sf/dipole/sf (gko) <pair_dipole>`
   * :doc:`lj/smooth (gko) <pair_lj_smooth>`
   * :doc:`lj/smooth/linear (ko) <pair_lj_smooth_linear>`
   * :doc:`lj/switch3/coulgauss/long (ko) <pair_lj_switch3_coulgauss_long>`
   * :doc:`lj96/cut (gko) <pair_lj96>`
   * :doc:`local/density (o) <pair_local_density>`
   * :doc:`lubricate (o) <pair_lubricate>`
   * :doc:`lubricate/poly (o) <pair_lubricate>`
   * :doc:`lubricateU <pair_lubricateU>`
   * :doc:`lubricateU/poly <pair_lubricateU>`
   * :doc:`mbx <pair_mbx>`
   * :doc:`mdpd (go) <pair_mesodpd>`
   * :doc:`mdpd/rhosum (o) <pair_mesodpd>`
   * :doc:`meam (k) <pair_meam>`
   * :doc:`meam/ms (k) <pair_meam>`
   * :doc:`meam/spline (o) <pair_meam_spline>`
   * :doc:`meam/sw/spline (o) <pair_meam_sw_spline>`
   * :doc:`mesocnt <pair_mesocnt>`
   * :doc:`mesocnt/viscous <pair_mesocnt>`
   * :doc:`mesomem/dipole (ko) <pair_mesomem_dipole>`
   * :doc:`mgpt <pair_mgpt>`
   * :doc:`mie/cut (gko) <pair_mie>`
   * :doc:`mliap (k) <pair_mliap>`
   * :doc:`mm3/switch3/coulgauss/long (ko) <pair_lj_switch3_coulgauss_long>`
   * :doc:`momb (ko) <pair_momb>`
   * :doc:`morse (gkot) <pair_morse>`
   * :doc:`morse/smooth/linear (ko) <pair_morse>`
   * :doc:`morse/soft (gko) <pair_fep_soft>`
   * :doc:`multi/lucy (o) <pair_multi_lucy>`
   * :doc:`multi/lucy/rx (k) <pair_multi_lucy_rx>`
   * :doc:`nb3b/harmonic (ko) <pair_nb3b>`
   * :doc:`nb3b/screened <pair_nb3b>`
   * :doc:`nm/cut (gko) <pair_nm>`
   * :doc:`nm/cut/coul/cut (gko) <pair_nm>`
   * :doc:`nm/cut/coul/long (gko) <pair_nm>`
   * :doc:`nm/cut/split (gko) <pair_nm>`
   * :doc:`oxdna/coaxstk <pair_oxdna>`
   * :doc:`oxdna/excv <pair_oxdna>`
   * :doc:`oxdna/hbond <pair_oxdna>`
   * :doc:`oxdna/stk <pair_oxdna>`
   * :doc:`oxdna/xstk <pair_oxdna>`
   * :doc:`oxdna2/coaxstk <pair_oxdna2>`
   * :doc:`oxdna2/dh <pair_oxdna2>`
   * :doc:`oxdna2/excv <pair_oxdna2>`
   * :doc:`oxdna2/hbond <pair_oxdna2>`
   * :doc:`oxdna2/stk <pair_oxdna2>`
   * :doc:`oxdna2/xstk <pair_oxdna2>`
   * :doc:`oxdna3/coaxstk <pair_oxdna3>`
   * :doc:`oxdna3/dh <pair_oxdna3>`
   * :doc:`oxdna3/excv <pair_oxdna3>`
   * :doc:`oxdna3/hbond <pair_oxdna3>`
   * :doc:`oxdna3/stk <pair_oxdna3>`
   * :doc:`oxdna3/xstk <pair_oxdna3>`
   * :doc:`oxrna2/excv <pair_oxrna2>`
   * :doc:`oxrna2/hbond <pair_oxrna2>`
   * :doc:`oxrna2/dh <pair_oxrna2>`
   * :doc:`oxrna2/stk <pair_oxrna2>`
   * :doc:`oxrna2/xstk <pair_oxrna2>`
   * :doc:`oxrna2/coaxstk <pair_oxrna2>`
   * :doc:`pace (ko) <pair_pace>`
   * :doc:`pace/extrapolation (k) <pair_pace>`
   * :doc:`pace/apip <pair_pace_apip>`
   * :doc:`pace/fast/apip <pair_pace_apip>`
   * :doc:`pace/precise/apip <pair_pace_apip>`
   * :doc:`pedone (ko) <pair_pedone>`
   * :doc:`pod (k) <pair_pod>`
   * :doc:`peri/eps <pair_peri>`
   * :doc:`peri/lps (o) <pair_peri>`
   * :doc:`peri/pmb (o) <pair_peri>`
   * :doc:`peri/ves <pair_peri>`
   * :doc:`polymorphic <pair_polymorphic>`
   * :doc:`quip <pair_quip>`
   * :doc:`reaxff (ko) <pair_reaxff>`
   * :doc:`rebo (io) <pair_airebo>`
   * :doc:`rebomos (o) <pair_rebomos>`
   * :doc:`resquared (go) <pair_resquared>`
   * :doc:`rheo <pair_rheo>`
   * :doc:`rheo/solid (o) <pair_rheo_solid>`
   * :doc:`runner <pair_runner>`
   * :doc:`saip/metal (t) <pair_saip_metal>`
   * :doc:`saip/metal/tmd (t) <pair_saip_metal>`
   * :doc:`sdpd/taitwater/isothermal (o) <pair_sdpd_taitwater_isothermal>`
   * :doc:`smatb (o) <pair_smatb>`
   * :doc:`smatb/single (o) <pair_smatb>`
   * :doc:`smd/hertz <pair_smd_hertz>`
   * :doc:`smd/tlsph <pair_smd_tlsph>`
   * :doc:`smd/tri_surface <pair_smd_triangulated_surface>`
   * :doc:`smd/ulsph <pair_smd_ulsph>`
   * :doc:`smtbq <pair_smtbq>`
   * :doc:`snap (ik) <pair_snap>`
   * :doc:`soft (gko) <pair_soft>`
   * :doc:`sph/heatconduction (go) <pair_sph_heatconduction>`
   * :doc:`sph/idealgas (o) <pair_sph_idealgas>`
   * :doc:`sph/lj (go) <pair_sph_lj>`
   * :doc:`sph/rhosum (o) <pair_sph_rhosum>`
   * :doc:`sph/taitwater (go) <pair_sph_taitwater>`
   * :doc:`sph/taitwater/morris (o) <pair_sph_taitwater_morris>`
   * :doc:`spin/dipole/cut <pair_spin_dipole>`
   * :doc:`spin/dipole/long <pair_spin_dipole>`
   * :doc:`spin/dmi <pair_spin_dmi>`
   * :doc:`spin/exchange <pair_spin_exchange>`
   * :doc:`spin/exchange/biquadratic <pair_spin_exchange>`
   * :doc:`spin/magelec <pair_spin_magelec>`
   * :doc:`spin/neel <pair_spin_neel>`
   * :doc:`srp <pair_srp>`
   * :doc:`srp/react <pair_srp>`
   * :doc:`surf/granular <pair_surf_granular>`
   * :doc:`sw (giko) <pair_sw>`
   * :doc:`sw/angle/table (ko) <pair_sw_angle_table>`
   * :doc:`sw/mod (ko) <pair_sw>`
   * :doc:`table/rx (k) <pair_table_rx>`
   * :doc:`tdpd (o) <pair_mesodpd>`
   * :doc:`tersoff (giko) <pair_tersoff>`
   * :doc:`tersoff/mod (gko) <pair_tersoff_mod>`
   * :doc:`tersoff/mod/c (ko) <pair_tersoff_mod>`
   * :doc:`tersoff/table (o) <pair_tersoff>`
   * :doc:`tersoff/zbl (gko) <pair_tersoff_zbl>`
   * :doc:`thole (o) <pair_thole>`
   * :doc:`threebody/table (ko) <pair_threebody_table>`
   * :doc:`tip4p/cut (ko) <pair_coul>`
   * :doc:`tip4p/long (ko) <pair_coul>`
   * :doc:`tip4p/long/soft (ko) <pair_fep_soft>`
   * :doc:`tri/lj (o) <pair_tri_lj>`
   * :doc:`ufm (gkot) <pair_ufm>`
   * :doc:`uf3 (ko) <pair_uf3>`
   * :doc:`vashishta (gko) <pair_vashishta>`
   * :doc:`vashishta/table (ko) <pair_vashishta>`
   * :doc:`wf/cut (ko) <pair_wf_cut>`
   * :doc:`ylz (ko) <pair_ylz>`
   * :doc:`yukawa (gko) <pair_yukawa>`
   * :doc:`yukawa/colloid (gko) <pair_yukawa_colloid>`
   * :doc:`zbl (gko) <pair_zbl>`
