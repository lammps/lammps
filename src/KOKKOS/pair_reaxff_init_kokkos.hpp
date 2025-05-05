// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_REAXFF_INIT_KOKKOS_HPP
#define LMP_PAIR_REAXFF_INIT_KOKKOS_HPP

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairReaxFFKokkos<DeviceType>::PairReaxFFKokkos(LAMMPS *lmp) : PairReaxFF(lmp)
{
  respa_enable = 0;

  cut_nbsq = cut_hbsq = cut_bosq = 0.0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | Q_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_resize_bo = DAT::tdual_int_scalar("pair:resize_bo");
  d_resize_bo = k_resize_bo.view<DeviceType>();

  k_resize_hb = DAT::tdual_int_scalar("pair:resize_hb");
  d_resize_hb = k_resize_hb.view<DeviceType>();

  nmax = 0;
  maxbo = 1;
  maxhb = 1;

  k_error_flag = DAT::tdual_int_scalar("pair:error_flag");
  k_nbuf_local = DAT::tdual_int_scalar("pair:nbuf_local");

  MemKK::realloc_kokkos(d_torsion_pack,"reaxff:torsion_pack",1,2);
  MemKK::realloc_kokkos(d_angular_pack,"reaxff:angular_pack",1,2);

  k_count_angular_torsion = DAT::tdual_int_1d("PairReaxFF::count_angular_torsion",2);
  d_count_angular_torsion = k_count_angular_torsion.template view<DeviceType>();

  if (execution_space == Host) list_blocking_flag = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairReaxFFKokkos<DeviceType>::~PairReaxFFKokkos()
{
  if (copymode) return;

  DeAllocate_System(api->system);

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  memoryKK->destroy_kokkos(k_tmpid,tmpid);
  tmpid = nullptr;
  memoryKK->destroy_kokkos(k_tmpbo,tmpbo);
  tmpbo = nullptr;

  deallocate_views_of_views();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::deallocate_views_of_views()
{

  // deallocate views of views in serial to prevent race conditions

  for (int i = 0; i < (int)k_LR.extent(0); i++) {
    for (int j = 0; j < (int)k_LR.extent(1); j++) {
      k_LR.h_view(i,j).d_vdW    = {};
      k_LR.h_view(i,j).d_CEvd   = {};
      k_LR.h_view(i,j).d_ele    = {};
      k_LR.h_view(i,j).d_CEclmb = {};
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::allocate()
{
  int n = atom->ntypes;

  k_params_sing = Kokkos::DualView<params_sing*,typename DeviceType::array_layout,DeviceType>
    ("PairReaxFF::params_sing",n+1);
  paramssing = k_params_sing.template view<DeviceType>();

  k_params_twbp = Kokkos::DualView<params_twbp**,typename DeviceType::array_layout,DeviceType>
    ("PairReaxFF::params_twbp",n+1,n+1);
  paramstwbp = k_params_twbp.template view<DeviceType>();

  k_params_thbp = Kokkos::DualView<params_thbp***,typename DeviceType::array_layout,DeviceType>
    ("PairReaxFF::params_thbp",n+1,n+1,n+1);
  paramsthbp = k_params_thbp.template view<DeviceType>();

  k_params_fbp = Kokkos::DualView<params_fbp****,typename DeviceType::array_layout,DeviceType>
    ("PairReaxFF::params_fbp",n+1,n+1,n+1,n+1);
  paramsfbp = k_params_fbp.template view<DeviceType>();

  k_params_hbp = Kokkos::DualView<params_hbp***,typename DeviceType::array_layout,DeviceType>
    ("PairReaxFF::params_hbp",n+1,n+1,n+1);
  paramshbp = k_params_hbp.template view<DeviceType>();

  k_tap = DAT::tdual_ffloat_1d("pair:tap",8);
  d_tap = k_tap.template view<DeviceType>();
  h_tap = k_tap.h_view;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::init_style()
{
  PairReaxFF::init_style();
  if (fix_reaxff) modify->delete_fix(fix_id); // not needed in the Kokkos version
  fix_reaxff = nullptr;

  acks2_flag = api->system->acks2_flag;
  if (acks2_flag) {
    auto ifix = modify->get_fix_by_style("^acks2/reax").front();
    if (!ifix->kokkosable)
      error->all(FLERR,"Must use Kokkos version of acks2/reaxff with pair reaxff/kk");
    if (ifix->execution_space == Host) {
      auto k_s = ((FixACKS2ReaxFFLegacyKokkos<LMPHostType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    } else {
      auto k_s = ((FixACKS2ReaxFFLegacyKokkos<LMPDeviceType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    }
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list with pair style reaxff/kk");

  need_dup = lmp->kokkos->need_dup<DeviceType>();

  allocate();
  setup();
  init_md();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::setup()
{
  int i,j,k,m;
  int n = atom->ntypes;

  // general parameters
  for (i = 0; i < 39; i ++)
    gp[i] = api->system->reax_param.gp.l[i];

  p_boc1 = gp[0];
  p_boc2 = gp[1];

  // vdw parameters
  vdwflag = api->system->reax_param.gp.vdw_type;
  lgflag = api->control->lgflag;

  // atom, bond, angle, dihedral, H-bond specific parameters
  two_body_parameters *twbp;

  // valence angle (3-body) parameters
  three_body_header *thbh;
  three_body_parameters *thbp;

  // torsion angle (4-body) parameters
  four_body_header *fbh;
  four_body_parameters *fbp;

  // hydrogen bond parameters
  hbond_parameters *hbp;

  for (i = 1; i <= n; i++) {
    if (map[i] == -1) continue;

    // general
    k_params_sing.h_view(i).mass = api->system->reax_param.sbp[map[i]].mass;

    // polarization
    k_params_sing.h_view(i).chi = api->system->reax_param.sbp[map[i]].chi;
    k_params_sing.h_view(i).eta = api->system->reax_param.sbp[map[i]].eta;

    // bond order
    k_params_sing.h_view(i).r_s = api->system->reax_param.sbp[map[i]].r_s;
    k_params_sing.h_view(i).r_p = api->system->reax_param.sbp[map[i]].r_p;
    k_params_sing.h_view(i).r_pp = api->system->reax_param.sbp[map[i]].r_pp;
    k_params_sing.h_view(i).valency = api->system->reax_param.sbp[map[i]].valency;
    k_params_sing.h_view(i).valency_val = api->system->reax_param.sbp[map[i]].valency_val;
    k_params_sing.h_view(i).valency_boc = api->system->reax_param.sbp[map[i]].valency_boc;
    k_params_sing.h_view(i).valency_e = api->system->reax_param.sbp[map[i]].valency_e;
    k_params_sing.h_view(i).nlp_opt = api->system->reax_param.sbp[map[i]].nlp_opt;

    // multibody
    k_params_sing.h_view(i).p_lp2 = api->system->reax_param.sbp[map[i]].p_lp2;
    k_params_sing.h_view(i).p_ovun2 = api->system->reax_param.sbp[map[i]].p_ovun2;
    k_params_sing.h_view(i).p_ovun5 = api->system->reax_param.sbp[map[i]].p_ovun5;

    // angular
    k_params_sing.h_view(i).p_val3 = api->system->reax_param.sbp[map[i]].p_val3;
    k_params_sing.h_view(i).p_val5 = api->system->reax_param.sbp[map[i]].p_val5;

    // hydrogen bond
    k_params_sing.h_view(i).p_hbond = api->system->reax_param.sbp[map[i]].p_hbond;

    // acks2
    k_params_sing.h_view(i).bcut_acks2 = api->system->reax_param.sbp[map[i]].bcut_acks2;

    for (j = 1; j <= n; j++) {
      if (map[j] == -1) continue;

      twbp = &(api->system->reax_param.tbp[map[i]][map[j]]);

      // vdW
      k_params_twbp.h_view(i,j).gamma = twbp->gamma;
      k_params_twbp.h_view(i,j).gamma_w = twbp->gamma_w;
      k_params_twbp.h_view(i,j).alpha = twbp->alpha;
      k_params_twbp.h_view(i,j).r_vdw = twbp->r_vdw;
      k_params_twbp.h_view(i,j).epsilon = twbp->D;
      k_params_twbp.h_view(i,j).acore = twbp->acore;
      k_params_twbp.h_view(i,j).ecore = twbp->ecore;
      k_params_twbp.h_view(i,j).rcore = twbp->rcore;
      k_params_twbp.h_view(i,j).lgre = twbp->lgre;
      k_params_twbp.h_view(i,j).lgcij = twbp->lgcij;

      // bond order
      k_params_twbp.h_view(i,j).r_s = twbp->r_s;
      k_params_twbp.h_view(i,j).r_p = twbp->r_p;
      k_params_twbp.h_view(i,j).r_pp = twbp->r_pp;
      k_params_twbp.h_view(i,j).p_bo1 = twbp->p_bo1;
      k_params_twbp.h_view(i,j).p_bo2 = twbp->p_bo2;
      k_params_twbp.h_view(i,j).p_bo3 = twbp->p_bo3;
      k_params_twbp.h_view(i,j).p_bo4 = twbp->p_bo4;
      k_params_twbp.h_view(i,j).p_bo5 = twbp->p_bo5;
      k_params_twbp.h_view(i,j).p_bo6 = twbp->p_bo6;
      k_params_twbp.h_view(i,j).p_boc3 = twbp->p_boc3;
      k_params_twbp.h_view(i,j).p_boc4 = twbp->p_boc4;
      k_params_twbp.h_view(i,j).p_boc5 = twbp->p_boc5;
      k_params_twbp.h_view(i,j).ovc = twbp->ovc;
      k_params_twbp.h_view(i,j).v13cor = twbp->v13cor;

      // bond energy
      k_params_twbp.h_view(i,j).p_be1 = twbp->p_be1;
      k_params_twbp.h_view(i,j).p_be2 = twbp->p_be2;
      k_params_twbp.h_view(i,j).De_s = twbp->De_s;
      k_params_twbp.h_view(i,j).De_p = twbp->De_p;
      k_params_twbp.h_view(i,j).De_pp = twbp->De_pp;

      // multibody
      k_params_twbp.h_view(i,j).p_ovun1 = twbp->p_ovun1;

      for (k = 1; k <= n; k++) {
        if (map[k] == -1) continue;

        // Angular
        thbh = &(api->system->reax_param.thbp[map[i]][map[j]][map[k]]);
        thbp = &(thbh->prm[0]);
        k_params_thbp.h_view(i,j,k).cnt = thbh->cnt;
        k_params_thbp.h_view(i,j,k).theta_00 = thbp->theta_00;
        k_params_thbp.h_view(i,j,k).p_val1 = thbp->p_val1;
        k_params_thbp.h_view(i,j,k).p_val2 = thbp->p_val2;
        k_params_thbp.h_view(i,j,k).p_val4 = thbp->p_val4;
        k_params_thbp.h_view(i,j,k).p_val7 = thbp->p_val7;
        k_params_thbp.h_view(i,j,k).p_pen1 = thbp->p_pen1;
        k_params_thbp.h_view(i,j,k).p_coa1 = thbp->p_coa1;

        // Hydrogen Bond
        hbp = &(api->system->reax_param.hbp[map[i]][map[j]][map[k]]);
        k_params_hbp.h_view(i,j,k).p_hb1 = hbp->p_hb1;
        k_params_hbp.h_view(i,j,k).p_hb2 = hbp->p_hb2;
        k_params_hbp.h_view(i,j,k).p_hb3 = hbp->p_hb3;
        k_params_hbp.h_view(i,j,k).r0_hb = hbp->r0_hb;

        for (m = 1; m <= n; m++) {
          if (map[m] == -1) continue;

          // Torsion
          fbh = &(api->system->reax_param.fbp[map[i]][map[j]][map[k]][map[m]]);
          fbp = &(fbh->prm[0]);
          k_params_fbp.h_view(i,j,k,m).p_tor1 = fbp->p_tor1;
          k_params_fbp.h_view(i,j,k,m).p_cot1 = fbp->p_cot1;
          k_params_fbp.h_view(i,j,k,m).V1 = fbp->V1;
          k_params_fbp.h_view(i,j,k,m).V2 = fbp->V2;
          k_params_fbp.h_view(i,j,k,m).V3 = fbp->V3;
        }
      }
    }
  }
  k_params_sing.template modify<LMPHostType>();
  k_params_twbp.template modify<LMPHostType>();
  k_params_thbp.template modify<LMPHostType>();
  k_params_fbp.template modify<LMPHostType>();
  k_params_hbp.template modify<LMPHostType>();

  // cutoffs
  cut_nbsq = api->control->nonb_cut * api->control->nonb_cut;
  cut_hbsq = api->control->hbond_cut * api->control->hbond_cut;
  cut_bosq = api->control->bond_cut * api->control->bond_cut;

  // bond order cutoffs
  bo_cut = 0.01 * gp[29];
  thb_cut = api->control->thb_cut;
  thb_cutsq = 0.000010; //thb_cut*thb_cut;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    allocate_array();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::init_md()
{
  // init_taper()
  F_FLOAT d1, d7, swa, swa2, swa3, swb, swb2, swb3;
  LR_lookup_table ** & LR = api->system->LR;

  swa = api->control->nonb_low;
  swb = api->control->nonb_cut;
  enobondsflag = api->control->enobondsflag;

  if ((fabs(swa) > 0.01) && (comm->me == 0))
    error->warning(FLERR, "Non-zero lower Taper-radius cutoff");

  if (swb < 0.0) {
    error->all(FLERR,"Negative upper Taper-radius cutoff");
  } else if ((swb < 5.0) && (comm->me ==0))
    error->warning(FLERR,"Very low Taper-radius cutoff: {}\n", swb);

  d1 = swb - swa;
  d7 = powint(d1,7);
  swa2 = swa * swa;
  swa3 = swa * swa2;
  swb2 = swb * swb;
  swb3 = swb * swb2;

  k_tap.h_view(7) = 20.0/d7;
  k_tap.h_view(6) = -70.0 * (swa + swb) / d7;
  k_tap.h_view(5) =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  k_tap.h_view(4) = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3) / d7;
  k_tap.h_view(3) = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3) / d7;
  k_tap.h_view(2) =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  k_tap.h_view(1) = 140.0 * swa3 * swb3 / d7;
  k_tap.h_view(0) = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
                     7.0*swa*swb3*swb3 + swb3*swb3*swb) / d7;

  k_tap.template modify<LMPHostType>();
  k_tap.template sync<DeviceType>();


  if (api->control->tabulate) {
    int ntypes = atom->ntypes;

    Init_Lookup_Tables();
    deallocate_views_of_views();
    k_LR = tdual_LR_lookup_table_kk_2d("lookup:LR",ntypes+1,ntypes+1);

    for (int i = 1; i <= ntypes; ++i) {
      if (map[i] == -1) continue;
      for (int j = i; j <= ntypes; ++j) {
        if (map[j] == -1) continue;
        int n = LR[i][j].n;
        if (n == 0) continue;
        k_LR.h_view(i,j).dx     = LR[i][j].dx;
        k_LR.h_view(i,j).inv_dx = LR[i][j].inv_dx;

        typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d k_vdW    = typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d("lookup:LR[i,j].vdW",n);
        typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d k_CEvd   = typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d("lookup:LR[i,j].CEvd",n);
        typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d k_ele    = typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d("lookup:LR[i,j].ele",n);
        typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d k_CEclmb = typename LR_lookup_table_kk<DeviceType>::tdual_cubic_spline_coef_1d("lookup:LR[i,j].CEclmb",n);

        k_LR.h_view(i,j).d_vdW    = k_vdW.template view<DeviceType>();
        k_LR.h_view(i,j).d_CEvd   = k_CEvd.template view<DeviceType>();
        k_LR.h_view(i,j).d_ele    = k_ele.template view<DeviceType>();
        k_LR.h_view(i,j).d_CEclmb = k_CEclmb.template view<DeviceType>();

        for (int k = 0; k < n; k++) {
          k_vdW.h_view(k)    = LR[i][j].vdW[k];
          k_CEvd.h_view(k)   = LR[i][j].CEvd[k];
          k_ele.h_view(k)    = LR[i][j].ele[k];
          k_CEclmb.h_view(k) = LR[i][j].CEclmb[k];
        }

        k_vdW.template modify<LMPHostType>();
        k_CEvd.template modify<LMPHostType>();
        k_ele.template modify<LMPHostType>();
        k_CEclmb.template modify<LMPHostType>();

        k_vdW.template sync<DeviceType>();
        k_CEvd.template sync<DeviceType>();
        k_ele.template sync<DeviceType>();
        k_CEclmb.template sync<DeviceType>();
      }
    }
    k_LR.template modify<LMPHostType>();
    k_LR.template sync<DeviceType>();

    Deallocate_Lookup_Tables();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairReaxFFKokkos<DeviceType>::Init_Lookup_Tables()
{
  int i, j, r;
  int num_atom_types;
  double dr;
  double *h, *fh, *fvdw, *fele, *fCEvd, *fCEclmb;
  double v0_vdw, v0_ele, vlast_vdw, vlast_ele;
  LR_lookup_table ** & LR = api->system->LR;

  /* initializations */
  v0_vdw = 0;
  v0_ele = 0;
  vlast_vdw = 0;
  vlast_ele = 0;

  num_atom_types = atom->ntypes;
  dr = api->control->nonb_cut / api->control->tabulate;
  h = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:h");
  fh = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:fh");
  fvdw = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:fvdw");
  fCEvd = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:fCEvd");
  fele = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:fele");
  fCEclmb = (double*)
    smalloc(api->control->error_ptr, (api->control->tabulate+2) * sizeof(double), "lookup:fCEclmb");

  LR = (LR_lookup_table**)
    scalloc(api->control->error_ptr, num_atom_types+1, sizeof(LR_lookup_table*), "lookup:LR");
  for (i = 0; i < num_atom_types+1; ++i)
    LR[i] = (LR_lookup_table*)
      scalloc(api->control->error_ptr, num_atom_types+1, sizeof(LR_lookup_table), "lookup:LR[i]");

  for (i = 1; i <= num_atom_types; ++i) {
    for (j = i; j <= num_atom_types; ++j) {
      LR[i][j].xmin = 0;
      LR[i][j].xmax = api->control->nonb_cut;
      LR[i][j].n = api->control->tabulate + 2;
      LR[i][j].dx = dr;
      LR[i][j].inv_dx = api->control->tabulate / api->control->nonb_cut;
      LR[i][j].y = (LR_data*)
        smalloc(api->control->error_ptr, LR[i][j].n * sizeof(LR_data), "lookup:LR[i,j].y");
      LR[i][j].H = (cubic_spline_coef*)
        smalloc(api->control->error_ptr, LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].H");
      LR[i][j].vdW = (cubic_spline_coef*)
        smalloc(api->control->error_ptr, LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].vdW");
      LR[i][j].CEvd = (cubic_spline_coef*)
        smalloc(api->control->error_ptr, LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].CEvd");
      LR[i][j].ele = (cubic_spline_coef*)
        smalloc(api->control->error_ptr, LR[i][j].n*sizeof(cubic_spline_coef),"lookup:LR[i,j].ele");
      LR[i][j].CEclmb = (cubic_spline_coef*)
        smalloc(api->control->error_ptr, LR[i][j].n*sizeof(cubic_spline_coef),
                 "lookup:LR[i,j].CEclmb");

      for (r = 1; r <= api->control->tabulate; ++r) {
        LR_vdW_Coulomb(i, j, r * dr, &(LR[i][j].y[r]));
        h[r] = LR[i][j].dx;
        fh[r] = LR[i][j].y[r].H;
        fvdw[r] = LR[i][j].y[r].e_vdW;
        fCEvd[r] = LR[i][j].y[r].CEvd;
        fele[r] = LR[i][j].y[r].e_ele;
        fCEclmb[r] = LR[i][j].y[r].CEclmb;
      }

      // init the start-end points
      h[r] = LR[i][j].dx;
      v0_vdw = LR[i][j].y[1].CEvd;
      v0_ele = LR[i][j].y[1].CEclmb;
      fh[r] = fh[r-1];
      fvdw[r] = fvdw[r-1];
      fCEvd[r] = fCEvd[r-1];
      fele[r] = fele[r-1];
      fCEclmb[r] = fCEclmb[r-1];
      vlast_vdw = fCEvd[r-1];
      vlast_ele = fele[r-1];

      Natural_Cubic_Spline(api->control->error_ptr, &h[1], &fh[1],
                            &(LR[i][j].H[1]), api->control->tabulate+1);

      Complete_Cubic_Spline(api->control->error_ptr, &h[1], &fvdw[1], v0_vdw, vlast_vdw,
                             &(LR[i][j].vdW[1]), api->control->tabulate+1);

      Natural_Cubic_Spline(api->control->error_ptr, &h[1], &fCEvd[1],
                            &(LR[i][j].CEvd[1]), api->control->tabulate+1);

      Complete_Cubic_Spline(api->control->error_ptr, &h[1], &fele[1], v0_ele, vlast_ele,
                             &(LR[i][j].ele[1]), api->control->tabulate+1);

      Natural_Cubic_Spline(api->control->error_ptr, &h[1], &fCEclmb[1],
                            &(LR[i][j].CEclmb[1]), api->control->tabulate+1);
    }
  }
  free(h);
  free(fh);
  free(fvdw);
  free(fCEvd);
  free(fele);
  free(fCEclmb);

  return 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::Deallocate_Lookup_Tables()
{
  int i, j;
  int ntypes;
  LR_lookup_table ** & LR = api->system->LR;

  ntypes = atom->ntypes;

  for (i = 0; i <= ntypes; ++i) {
    if (map[i] == -1) continue;
    for (j = i; j <= ntypes; ++j) {
      if (map[j] == -1) continue;
      if (LR[i][j].n) {
        sfree(LR[i][j].y);
        sfree(LR[i][j].H);
        sfree(LR[i][j].vdW);
        sfree(LR[i][j].CEvd);
        sfree(LR[i][j].ele);
        sfree(LR[i][j].CEclmb);
      }
    }
    sfree(LR[i]);
  }
  sfree(LR);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::LR_vdW_Coulomb(int i, int j, double r_ij, LR_data *lr)
{
  double p_vdW1 = api->system->reax_param.gp.l[28];
  double p_vdW1i = 1.0 / p_vdW1;
  double powr_vdW1, powgi_vdW1;
  double tmp, fn13, exp1, exp2;
  double Tap, dTap, dfn13;
  double dr3gamij_1, dr3gamij_3;
  double e_core, de_core;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  two_body_parameters *twbp;

  twbp = &(api->system->reax_param.tbp[map[i]][map[j]]);
  e_core = 0;
  de_core = 0;
  e_lg = de_lg = 0.0;

  /* calculate taper and its derivative */
  Tap = k_tap.h_view[7] * r_ij + k_tap.h_view[6];
  Tap = Tap * r_ij + k_tap.h_view[5];
  Tap = Tap * r_ij + k_tap.h_view[4];
  Tap = Tap * r_ij + k_tap.h_view[3];
  Tap = Tap * r_ij + k_tap.h_view[2];
  Tap = Tap * r_ij + k_tap.h_view[1];
  Tap = Tap * r_ij + k_tap.h_view[0];

  dTap = 7*k_tap.h_view[7] * r_ij + 6*k_tap.h_view[6];
  dTap = dTap * r_ij + 5*k_tap.h_view[5];
  dTap = dTap * r_ij + 4*k_tap.h_view[4];
  dTap = dTap * r_ij + 3*k_tap.h_view[3];
  dTap = dTap * r_ij + 2*k_tap.h_view[2];
  dTap += k_tap.h_view[1]/r_ij;

  /*vdWaals Calculations*/
  if (api->system->reax_param.gp.vdw_type==1 || api->system->reax_param.gp.vdw_type==3)
    { // shielding
      powr_vdW1 = pow(r_ij, p_vdW1);
      powgi_vdW1 = pow(1.0 / twbp->gamma_w, p_vdW1);

      fn13 = pow(powr_vdW1 + powgi_vdW1, p_vdW1i);
      exp1 = exp(twbp->alpha * (1.0 - fn13 / twbp->r_vdw));
      exp2 = exp(0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdw));

      lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);

      dfn13 = pow(powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) * pow(r_ij, p_vdW1-2.0);

      lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
        Tap * twbp->D * (twbp->alpha / twbp->r_vdw) * (exp1 - exp2) * dfn13;
    }
  else { // no shielding
    exp1 = exp(twbp->alpha * (1.0 - r_ij / twbp->r_vdw));
    exp2 = exp(0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdw));

    lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);
    lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
      Tap * twbp->D * (twbp->alpha / twbp->r_vdw) * (exp1 - exp2) / r_ij;
  }

  if (api->system->reax_param.gp.vdw_type==2 || api->system->reax_param.gp.vdw_type==3)
    { // inner wall
      e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij/twbp->rcore)));
      lr->e_vdW += Tap * e_core;

      de_core = -(twbp->acore/twbp->rcore) * e_core;
      lr->CEvd += dTap * e_core + Tap * de_core / r_ij;

      //  lg correction, only if lgvdw is yes
      if (api->control->lgflag) {
        r_ij5 = powint(r_ij, 5);
        r_ij6 = powint(r_ij, 6);
        re6 = powint(twbp->lgre, 6);
        e_lg = -(twbp->lgcij/(r_ij6 + re6));
        lr->e_vdW += Tap * e_lg;

        de_lg = -6.0 * e_lg *  r_ij5 / (r_ij6 + re6) ;
        lr->CEvd += dTap * e_lg + Tap * de_lg/r_ij;
      }
    }

  /* Coulomb calculations */
  dr3gamij_1 = (r_ij * r_ij * r_ij + twbp->gamma);
  dr3gamij_3 = cbrt(dr3gamij_1);
  tmp = Tap / dr3gamij_3;
  lr->H = EV_to_KCALpMOL * tmp;
  lr->e_ele = C_ele * tmp;

  lr->CEclmb = C_ele * (dTap -  Tap * r_ij / dr3gamij_1) / dr3gamij_3;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::allocate_array()
{
  // free scatterview memory
  if (need_dup) {
    dup_dDeltap_self = {};
    dup_total_bo     = {};
    dup_CdDelta      = {};
  } else {
    ndup_dDeltap_self = {};
    ndup_total_bo     = {};
    ndup_CdDelta      = {};
  }

  if (cut_hbsq > 0.0) {
    MemKK::realloc_kokkos(d_hb_num,"reaxff/kk:hb_num",nmax);
    MemKK::realloc_kokkos(d_hb_list,"reaxff/kk:hb_list", nmax, maxhb);
  }
  MemKK::realloc_kokkos(d_bo_num,"reaxff/kk:bo_num",nmax);
  MemKK::realloc_kokkos(d_bo_list,"reaxff/kk:bo_list", nmax, maxbo);

  MemKK::realloc_kokkos(d_BO,"reaxff/kk:BO",nmax,maxbo);
  MemKK::realloc_kokkos(d_BO_s,"reaxff/kk:BO",nmax,maxbo);
  MemKK::realloc_kokkos(d_BO_pi,"reaxff/kk:BO_pi",nmax,maxbo);
  MemKK::realloc_kokkos(d_BO_pi2,"reaxff/kk:BO_pi2",nmax,maxbo);

  MemKK::realloc_kokkos(d_dln_BOp_pi,"reaxff/kk:d_dln_BOp_pi",nmax,maxbo);
  MemKK::realloc_kokkos(d_dln_BOp_pi2,"reaxff/kk:d_dln_BOp_pi2",nmax,maxbo);

  MemKK::realloc_kokkos(d_C1dbo,"reaxff/kk:d_C1dbo",nmax,maxbo);
  MemKK::realloc_kokkos(d_C2dbo,"reaxff/kk:d_C2dbo",nmax,maxbo);
  MemKK::realloc_kokkos(d_C3dbo,"reaxff/kk:d_C3dbo",nmax,maxbo);

  MemKK::realloc_kokkos(d_C1dbopi,"reaxff/kk:d_C1dbopi",nmax,maxbo);
  MemKK::realloc_kokkos(d_C2dbopi,"reaxff/kk:d_C2dbopi",nmax,maxbo);
  MemKK::realloc_kokkos(d_C3dbopi,"reaxff/kk:d_C3dbopi",nmax,maxbo);
  MemKK::realloc_kokkos(d_C4dbopi,"reaxff/kk:d_C4dbopi",nmax,maxbo);

  MemKK::realloc_kokkos(d_C1dbopi2,"reaxff/kk:d_C1dbopi2",nmax,maxbo);
  MemKK::realloc_kokkos(d_C2dbopi2,"reaxff/kk:d_C2dbopi2",nmax,maxbo);
  MemKK::realloc_kokkos(d_C3dbopi2,"reaxff/kk:d_C3dbopi2",nmax,maxbo);
  MemKK::realloc_kokkos(d_C4dbopi2,"reaxff/kk:d_C4dbopi2",nmax,maxbo);

  MemKK::realloc_kokkos(d_dBOp,"reaxff/kk:dBOp",nmax,maxbo);

  MemKK::realloc_kokkos(d_dDeltap_self,"reaxff/kk:dDeltap_self",nmax,3);
  MemKK::realloc_kokkos(d_Deltap_boc,"reaxff/kk:Deltap_boc",nmax);
  MemKK::realloc_kokkos(d_Deltap,"reaxff/kk:Deltap",nmax);
  MemKK::realloc_kokkos(d_total_bo,"reaxff/kk:total_bo",nmax);

  MemKK::realloc_kokkos(d_Cdbo,"reaxff/kk:Cdbo",nmax,maxbo);
  MemKK::realloc_kokkos(d_Cdbopi,"reaxff/kk:Cdbopi",nmax,maxbo);
  MemKK::realloc_kokkos(d_Cdbopi2,"reaxff/kk:Cdbopi2",nmax,maxbo);

  MemKK::realloc_kokkos(d_Delta,"reaxff/kk:Delta",nmax);
  MemKK::realloc_kokkos(d_Delta_boc,"reaxff/kk:Delta_boc",nmax);
  MemKK::realloc_kokkos(d_dDelta_lp,"reaxff/kk:dDelta_lp",nmax);
  MemKK::realloc_kokkos(d_Delta_lp,"reaxff/kk:Delta_lp",nmax);
  MemKK::realloc_kokkos(d_Delta_lp_temp,"reaxff/kk:Delta_lp_temp",nmax);
  MemKK::realloc_kokkos(d_CdDelta,"reaxff/kk:CdDelta",nmax);
  MemKK::realloc_kokkos(d_sum_ovun,"reaxff/kk:sum_ovun",nmax,3);

  // FixReaxFFBonds
  MemKK::realloc_kokkos(d_abo,"reaxff/kk:abo",nmax,maxbo);
  MemKK::realloc_kokkos(d_neighid,"reaxff/kk:neighid",nmax,maxbo);
  MemKK::realloc_kokkos(d_numneigh_bonds,"reaxff/kk:numneigh_bonds",nmax);

  // ComputeAngular intermediates
  MemKK::realloc_kokkos(d_angular_intermediates,"reaxff/kk:angular_intermediates",nmax,4);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxZero, const int &n) const {
  d_total_bo(n) = 0.0;
  d_CdDelta(n) = 0.0;
  d_bo_num(n) = 0.0;
  d_hb_num(n) = 0.0;
  for (int j = 0; j < 3; j++)
    d_dDeltap_self(n,j) = 0.0;
}


#endif
