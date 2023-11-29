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

/* ----------------------------------------------------------------------
   Contributing authors: Ray Shan (SNL), Stan Moore (SNL),
     Evan Weinberg (NVIDIA)

   Nicholas Curtis (AMD), Leopold Grinberd (AMD), and Gina Sitaraman (AMD):
     - Reduced math overhead: enabled specialized calls (e.g., cbrt for a
         cube root instead of pow) and use power/exponential laws to reduce the
         number of exponentials evaluated, etc.
     - Added blocking to the Torsion and (optionally) BuildLists kernels, to
         reduce thread divergence on GPUs
     - Added preview to BuildLists kernels along with full version
------------------------------------------------------------------------- */

#include "pair_reaxff_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "fix_acks2_reaxff_kokkos.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_special.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "reaxff_api.h"
#include "reaxff_defs.h"

#include <cmath>

/* ---------------------------------------------------------------------- */

using namespace ReaxFF;

namespace LAMMPS_NS {
using namespace MathConst;
using namespace MathSpecial;

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

  // deallocate views of views in serial to prevent race condition in profiling tools

  for (int i = 0; i < (int)k_LR.extent(0); i++) {
    for (int j = 0; j < (int)k_LR.extent(1); j++) {
      k_LR.h_view(i,j).d_vdW    = decltype(k_LR.h_view(i,j).d_vdW   )();
      k_LR.h_view(i,j).d_CEvd   = decltype(k_LR.h_view(i,j).d_CEvd  )();
      k_LR.h_view(i,j).d_ele    = decltype(k_LR.h_view(i,j).d_ele   )();
      k_LR.h_view(i,j).d_CEclmb = decltype(k_LR.h_view(i,j).d_CEclmb)();
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
      auto k_s = ((FixACKS2ReaxFFKokkos<LMPHostType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    } else {
      auto k_s = ((FixACKS2ReaxFFKokkos<LMPDeviceType>*) ifix)->get_s();
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
    k_params_sing.h_view(i).r_pi = api->system->reax_param.sbp[map[i]].r_pi;
    k_params_sing.h_view(i).r_pi2 = api->system->reax_param.sbp[map[i]].r_pi_pi;
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
      k_params_twbp.h_view(i,j).r_vdw = twbp->r_vdW;
      k_params_twbp.h_view(i,j).epsilon = twbp->D;
      k_params_twbp.h_view(i,j).acore = twbp->acore;
      k_params_twbp.h_view(i,j).ecore = twbp->ecore;
      k_params_twbp.h_view(i,j).rcore = twbp->rcore;
      k_params_twbp.h_view(i,j).lgre = twbp->lgre;
      k_params_twbp.h_view(i,j).lgcij = twbp->lgcij;

      // bond order
      k_params_twbp.h_view(i,j).r_s = twbp->r_s;
      k_params_twbp.h_view(i,j).r_pi = twbp->r_p;
      k_params_twbp.h_view(i,j).r_pi2 = twbp->r_pp;
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

  if (fabs(swa) > 0.01)
    error->warning(FLERR,"Warning: non-zero lower Taper-radius cutoff");

  if (swb < 0)
    error->one(FLERR,"Negative upper Taper-radius cutoff");
  else if (swb < 5)
    error->one(FLERR,"Warning: very low Taper-radius cutoff: {}\n", swb);

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
    k_LR = tdual_LR_lookup_table_kk_2d("lookup:LR",ntypes+1,ntypes+1);
    d_LR = k_LR.template view<DeviceType>();

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
      if (map[i] == -1) continue;
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
      exp1 = exp(twbp->alpha * (1.0 - fn13 / twbp->r_vdW));
      exp2 = exp(0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW));

      lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);

      dfn13 = pow(powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) * pow(r_ij, p_vdW1-2.0);

      lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
        Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
    }
  else { // no shielding
    exp1 = exp(twbp->alpha * (1.0 - r_ij / twbp->r_vdW));
    exp2 = exp(0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW));

    lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);
    lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
      Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) / r_ij;
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
void PairReaxFFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  bocnt = hbcnt = 0;

  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  atomKK->sync(execution_space,datamask_read);
  k_params_sing.template sync<DeviceType>();
  k_params_twbp.template sync<DeviceType>();
  k_params_thbp.template sync<DeviceType>();
  k_params_fbp.template sync<DeviceType>();
  k_params_hbp.template sync<DeviceType>();

  if (eflag_either || vflag_either) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atomKK->nlocal;
  newton_pair = force->newton_pair;

  nn = list->inum;
  NN = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  if (acks2_flag) {
    auto ifix = modify->get_fix_by_style("^acks2/reax").front();
    if (ifix->execution_space == Host) {
      auto k_s = ((FixACKS2ReaxFFKokkos<LMPHostType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    } else {
      auto k_s = ((FixACKS2ReaxFFKokkos<LMPDeviceType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    }
  }

  // allocate duplicated memory
  if (need_dup) {
    dup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  if (eflag_global) {
    for (int i = 0; i < 14; i++)
      pvector[i] = 0.0;
  }

  EV_FLOAT_REAX ev;
  EV_FLOAT_REAX ev_all;

  // Polarization (self)
  if (neighflag == HALF) {
    if (eflag_global)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALF>>(0,inum),*this,ev);
    else if (eflag_atom)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALF>>(0,inum),*this);
  } else { //if (neighflag == HALFTHREAD) {
    if (eflag_global)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALFTHREAD>>(0,inum),*this,ev);
    else if (eflag_atom)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALFTHREAD>>(0,inum),*this);
  }
  ev_all += ev;
  pvector[13] = ev.ecoul;

  // LJ + Coulomb
  if (api->control->tabulate) {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALF,0>>(0,inum),*this);
    } else if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALFTHREAD,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALFTHREAD,0>>(0,inum),*this);
    }
  } else {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALF,0>>(0,inum),*this);
    } else if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALFTHREAD,1>,Kokkos::LaunchBounds<256,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALFTHREAD,0>,Kokkos::LaunchBounds<256,1>>(0,inum),*this);
    }
  }
  ev_all += ev;
  pvector[10] = ev.evdwl;
  pvector[11] = ev.ecoul;


  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    allocate_array();
  }

  // allocate duplicated memory
  if (need_dup) {
    dup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_dDeltap_self);
    dup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_total_bo);
  } else {
    ndup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_dDeltap_self);
    ndup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_total_bo);
  }

  // Neighbor lists for bond and hbond

  // try, resize if necessary

  int resize = 1;
  while (resize) {
    resize = 0;

    k_resize_bo.h_view() = 0;
    k_resize_bo.modify<LMPHostType>();
    k_resize_bo.sync<DeviceType>();

    k_resize_hb.h_view() = 0;
    k_resize_hb.modify<LMPHostType>();
    k_resize_hb.sync<DeviceType>();

    // zero
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxZero>(0,nmax),*this);

    if (execution_space == Host) { // CPU
      if (neighflag == HALF)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlocking<HALF>>(0,ignum),*this);
      else if (neighflag == HALFTHREAD)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlocking<HALFTHREAD>>(0,ignum),*this);
    } else {
      if (list_blocking_flag) {
        if (neighflag == HALF)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlockingPreview<HALF>>(0,ignum),*this);
        else if (neighflag == HALFTHREAD)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlockingPreview<HALFTHREAD>>(0,ignum),*this);
      } else {
        if (neighflag == HALF)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfPreview<HALF>>(0,ignum),*this);
        else if (neighflag == HALFTHREAD)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfPreview<HALFTHREAD>>(0,ignum),*this);
      }
    }

    k_resize_bo.modify<DeviceType>();
    k_resize_bo.sync<LMPHostType>();
    int resize_bo = k_resize_bo.h_view();
    if (resize_bo) maxbo = MAX(maxbo+MAX(1,maxbo*0.1),resize_bo);

    k_resize_hb.modify<DeviceType>();
    k_resize_hb.sync<LMPHostType>();
    int resize_hb = k_resize_hb.h_view();
    if (resize_hb) maxhb = MAX(maxhb+MAX(1,maxhb*0.1),resize_hb);

    resize = resize_bo || resize_hb;
    if (resize) {
      allocate_array();
      if (need_dup) {
        dup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_dDeltap_self);
        dup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_total_bo);
      } else {
        ndup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_dDeltap_self);
        ndup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_total_bo);
      }
    }
  }

  if (execution_space != Host) // GPU
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsFull>(0,ignum),*this);

  // allocate duplicated memory
  if (need_dup)
    dup_CdDelta = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_CdDelta);
  else
    ndup_CdDelta = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_CdDelta);

  // reduction over duplicated memory
  if (need_dup)
    Kokkos::Experimental::contribute(d_total_bo, dup_total_bo); // needed in BondOrder1

  // Bond order
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder1>(0,ignum),*this);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder2>(0,ignum),*this);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder3>(0,ignum),*this);

  // Bond energy
  if (neighflag == HALF) {
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALF,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALF,0>>(0,inum),*this);
    ev_all += ev;
    pvector[0] = ev.evdwl;
  } else { //if (neighflag == HALFTHREAD) {
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALFTHREAD,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALFTHREAD,0>>(0,inum),*this);
    ev_all += ev;
    pvector[0] = ev.evdwl;
  }

  // Multi-body corrections
  if (neighflag == HALF) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti1>(0,inum),*this);
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALF,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALF,0>>(0,inum),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti1>(0,inum),*this);
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALFTHREAD,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALFTHREAD,0>>(0,inum),*this);
    ev_all += ev;
  }
  pvector[2] = ev.ereax[0];
  pvector[1] = ev.ereax[1]+ev.ereax[2];
  pvector[3] = 0.0;
  ev_all.evdwl += ev.ereax[0] + ev.ereax[1] + ev.ereax[2];

  int count_angular = 0;
  int count_torsion = 0;

  auto& h_count_angular_torsion = k_count_angular_torsion.h_view;
  h_count_angular_torsion(0) = 0;
  h_count_angular_torsion(1) = 0;
  k_count_angular_torsion.template modify<LMPHostType>();
  k_count_angular_torsion.template sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxCountAngularTorsion<false> >(0,inum),*this);

  k_count_angular_torsion.template modify<DeviceType>();
  k_count_angular_torsion.template sync<LMPHostType>();
  count_angular = h_count_angular_torsion(0);
  count_torsion = h_count_angular_torsion(1);

  if (count_angular > (int)d_angular_pack.extent(0))
    MemKK::realloc_kokkos(d_angular_pack,"reaxff:angular_pack",(int)(count_angular * 1.1),2);
  if (count_torsion > (int)d_torsion_pack.extent(0))
    MemKK::realloc_kokkos(d_torsion_pack,"reaxff:torsion_pack",(int)(count_torsion * 1.1),2);

  // need to zero to re-count
  h_count_angular_torsion(0) = 0;
  h_count_angular_torsion(1) = 0;
  k_count_angular_torsion.template modify<LMPHostType>();
  k_count_angular_torsion.template sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxCountAngularTorsion<true>>(0,inum),*this);

  // no need to re-sync count_angular, count_torsion

  // Angular
  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALF,1>>(0,count_angular),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALF,0>>(0,count_angular),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALFTHREAD,1>>(0,count_angular),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALFTHREAD,0>>(0,count_angular),*this);
    ev_all += ev;
  }
  pvector[4] = ev.ereax[3];
  pvector[5] = ev.ereax[4];
  pvector[6] = ev.ereax[5];
  ev_all.evdwl += ev.ereax[3] + ev.ereax[4] + ev.ereax[5];

  // Torsion
  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALF,1>>(0,count_torsion),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALF,0>>(0,count_torsion),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALFTHREAD,1>>(0,count_torsion),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALFTHREAD,0>>(0,count_torsion),*this);
    ev_all += ev;
  }
  pvector[8] = ev.ereax[6];
  pvector[9] = ev.ereax[7];
  ev_all.evdwl += ev.ereax[6] + ev.ereax[7];

  // Hydrogen Bond
  if (cut_hbsq > 0.0) {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALF,0>>(0,inum),*this);
      ev_all += ev;
    } else { //if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALFTHREAD,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALFTHREAD,0>>(0,inum),*this);
      ev_all += ev;
    }
  }
  pvector[7] = ev.ereax[8];
  ev_all.evdwl += ev.ereax[8];

  // reduction over duplicated memory
  if (need_dup) {
    Kokkos::Experimental::contribute(d_dDeltap_self, dup_dDeltap_self); // needed in ComputeBond2
    Kokkos::Experimental::contribute(d_CdDelta, dup_CdDelta); // needed in ComputeBond2
  }

  // Bond force
  if (neighflag == HALF) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxUpdateBond<HALF>>(0,ignum),*this);

    if (vflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALF,1>>(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALF,0>>(0,ignum),*this);
    ev_all += ev;
    pvector[0] += ev.evdwl;
  } else { //if (neighflag == HALFTHREAD) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxUpdateBond<HALFTHREAD>>(0,ignum),*this);

    if (vflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALFTHREAD,1>>(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALFTHREAD,0>>(0,ignum),*this);
    ev_all += ev;
    pvector[0] += ev.evdwl;
  }

  // reduction over duplicated memory
  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) {
    eng_vdwl += ev_all.evdwl;
    eng_coul += ev_all.ecoul;
  }

  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (fixspecies_flag)
    FindBondSpecies();

  copymode = 0;

  // free scatterview memory
  if (need_dup) {
    dup_f            = decltype(dup_f)();
    dup_eatom        = decltype(dup_eatom)();
    dup_vatom        = decltype(dup_vatom)();
    dup_dDeltap_self = decltype(dup_dDeltap_self)();
    dup_total_bo     = decltype(dup_total_bo)();
    dup_CdDelta      = decltype(dup_CdDelta)();
  } else {
    ndup_f            = decltype(ndup_f)();
    ndup_eatom        = decltype(ndup_eatom)();
    ndup_vatom        = decltype(ndup_vatom)();
    ndup_dDeltap_self = decltype(ndup_dDeltap_self)();
    ndup_total_bo     = decltype(ndup_total_bo)();
    ndup_CdDelta      = decltype(ndup_CdDelta)();
  }

  d_neighbors = typename AT::t_neighbors_2d();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT qi = q(i);
  const F_FLOAT chi = paramssing(itype).chi;
  const F_FLOAT eta = paramssing(itype).eta;

  F_FLOAT epol = KCALpMOL_to_EV*(chi*qi+(eta/2.0)*qi*qi);

  /* energy due to coupling with kinetic energy potential */
  if (acks2_flag)
    epol += KCALpMOL_to_EV*qi*d_s[NN + i];

  if (eflag_global) ev.ecoul += epol;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,epol,0.0,0.0,0.0,0.0);
  if (eflag_atom) this->template e_tally_single<NEIGHFLAG>(ev,i,epol);
}

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG>(TagPairReaxComputePolar<NEIGHFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT powr_vdw, powgi_vdw, fn13, dfn13, exp1, exp2, etmp;
  F_FLOAT evdwl, fvdwl;
  evdwl = fvdwl = 0.0;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const F_FLOAT qi = q(i);
  const int itype = type(i);
  const tagint itag = tag(i);
  const int jnum = d_numneigh[i];

  F_FLOAT fxtmp, fytmp, fztmp;
  fxtmp = fytmp = fztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const tagint jtag = tag(j);
    const F_FLOAT qj = q(j);

    // skip half of the interactions
    if (j >= nlocal) {
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x(j,2) < ztmp) continue;
        if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
        if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
      }
    }

    const X_FLOAT delx = x(j,0) - xtmp;
    const X_FLOAT dely = x(j,1) - ytmp;
    const X_FLOAT delz = x(j,2) - ztmp;
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq > cut_nbsq) continue;
    const F_FLOAT rij = sqrt(rsq);

    // LJ energy/force
    F_FLOAT Tap = d_tap[7] * rij + d_tap[6];
    Tap = Tap * rij + d_tap[5];
    Tap = Tap * rij + d_tap[4];
    Tap = Tap * rij + d_tap[3];
    Tap = Tap * rij + d_tap[2];
    Tap = Tap * rij + d_tap[1];
    Tap = Tap * rij + d_tap[0];

    F_FLOAT dTap = 7*d_tap[7] * rij + 6*d_tap[6];
    dTap = dTap * rij + 5*d_tap[5];
    dTap = dTap * rij + 4*d_tap[4];
    dTap = dTap * rij + 3*d_tap[3];
    dTap = dTap * rij + 2*d_tap[2];
    dTap += d_tap[1]/rij;

    const F_FLOAT gamma_w = paramstwbp(itype,jtype).gamma_w;
    const F_FLOAT alpha = paramstwbp(itype,jtype).alpha;
    const F_FLOAT r_vdw = paramstwbp(itype,jtype).r_vdw;
    const F_FLOAT epsilon = paramstwbp(itype,jtype).epsilon;

    // shielding
    if (vdwflag == 1 || vdwflag == 3) {
      F_FLOAT tmp_var;
      tmp_var = pow(rij,gp[28]-2.0);
      powr_vdw = tmp_var*rij*rij;
      powgi_vdw = pow(1.0/gamma_w,gp[28]);
      dfn13 = pow(powr_vdw+powgi_vdw,1.0/gp[28]-1.0);
      fn13  = dfn13*(powr_vdw+powgi_vdw);
      dfn13 = dfn13*tmp_var;

      exp2 = exp(0.5*alpha*(1.0-fn13/r_vdw));
      exp1 = exp2*exp2;
      etmp = epsilon*(exp1-2.0*exp2);
      evdwl = Tap*etmp;
      fvdwl = dTap*etmp-Tap*epsilon*(alpha/r_vdw)*(exp1-exp2)*dfn13;
    } else {
      exp2 = exp(0.5*alpha*(1.0-rij/r_vdw));
      exp1 = exp2*exp2;
      etmp = epsilon*(exp1-2.0*exp2);
      evdwl = Tap*etmp;
      fvdwl = dTap*etmp-Tap*epsilon*(alpha/r_vdw)*(exp1-exp2)*rij;
    }
    // inner wall
    if (vdwflag == 2 || vdwflag == 3) {
      const F_FLOAT ecore = paramstwbp(itype,jtype).ecore;
      const F_FLOAT acore = paramstwbp(itype,jtype).acore;
      const F_FLOAT rcore = paramstwbp(itype,jtype).rcore;
      const F_FLOAT e_core = ecore*exp(acore*(1.0-(rij/rcore)));
      const F_FLOAT de_core = -(acore/rcore)*e_core;
      evdwl += Tap*e_core;
      fvdwl += dTap*e_core+Tap*de_core/rij;

      if (lgflag) {
        const F_FLOAT lgre = paramstwbp(itype,jtype).lgre;
        const F_FLOAT lgcij = paramstwbp(itype,jtype).lgcij;
        const F_FLOAT rij5 = rsq*rsq*rij;
        const F_FLOAT rij6 = rij5*rij;
        const F_FLOAT re6 = lgre*lgre*lgre*lgre*lgre*lgre;
        const F_FLOAT elg = -lgcij/(rij6+re6);
        const F_FLOAT delg = -6.0*elg*rij5/(rij6+re6);
        evdwl += Tap*elg;
        fvdwl += dTap*elg+Tap*delg/rij;
      }
    }

    // Coulomb energy/force
    const F_FLOAT shld = paramstwbp(itype,jtype).gamma;
    const F_FLOAT denom1 = rij * rij * rij + shld;
    const F_FLOAT denom3 = cbrt(denom1);
    F_FLOAT ecoul = C_ele * qi*qj*Tap/denom3;
    F_FLOAT fcoul = C_ele * qi*qj*(dTap-Tap*rij/denom1)/denom3;

    /* contribution to energy and gradients (atoms and cell)
     * due to geometry-dependent terms in the ACKS2
     * kinetic energy */
    if (acks2_flag) {

      /* kinetic energy terms */
      double xcut = 0.5 * (paramssing(itype).bcut_acks2
                          + paramssing(jtype).bcut_acks2);

      if (rij <= xcut) {
        const F_FLOAT d = rij / xcut;
        const F_FLOAT bond_softness = gp[34] * pow( d, 3.0 )
                                    * pow( 1.0 - d, 6.0 );

        if (bond_softness > 0.0) {
          /* Coulombic energy contribution */
          const F_FLOAT effpot_diff = d_s[NN + i]
                                    - d_s[NN + j];
          const F_FLOAT e_ele = -0.5 * KCALpMOL_to_EV * bond_softness
                                     * SQR( effpot_diff );

          ecoul += e_ele;

          /* forces contribution */
          F_FLOAT d_bond_softness;
          d_bond_softness = gp[34]
                          * 3.0 / xcut * pow( d, 2.0 )
                          * pow( 1.0 - d, 5.0 ) * (1.0 - 3.0 * d);
          d_bond_softness = -0.5 * d_bond_softness
                          * SQR( effpot_diff );
          d_bond_softness = KCALpMOL_to_EV * d_bond_softness
                          / rij;

          fcoul += d_bond_softness;
        }
      }
    }

    const F_FLOAT ftotal = fvdwl + fcoul;
    fxtmp += delx*ftotal;
    a_f(j,0) -= delx*ftotal;
    fytmp += dely*ftotal;
    a_f(j,1) -= dely*ftotal;
    fztmp += delz*ftotal;
    a_f(j,2) -= delz*ftotal;

    if (EVFLAG) {
      if (eflag_global) ev.evdwl += evdwl;
      if (eflag_global) ev.ecoul += ecoul;

      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl+ecoul,-ftotal,delx,dely,delz);
    }
  }

  a_f(i,0) += fxtmp;
  a_f(i,1) += fytmp;
  a_f(i,2) += fztmp;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const F_FLOAT qi = q(i);
  const int itype = type(i);
  const tagint itag = tag(i);
  const int jnum = d_numneigh[i];

  F_FLOAT fxtmp, fytmp, fztmp;
  fxtmp = fytmp = fztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const tagint jtag = tag(j);
    const F_FLOAT qj = q(j);

    // skip half of the interactions
    if (j >= nlocal) {
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x(j,2) < ztmp) continue;
        if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
        if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
      }
    }

    const X_FLOAT delx = x(j,0) - xtmp;
    const X_FLOAT dely = x(j,1) - ytmp;
    const X_FLOAT delz = x(j,2) - ztmp;
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq > cut_nbsq) continue;
    const F_FLOAT rij = sqrt(rsq);

    const int tmin  = MIN(itype, jtype);
    const int tmax  = MAX(itype, jtype);
    const LR_lookup_table_kk<DeviceType>& t = d_LR(tmin,tmax);


    /* Cubic Spline Interpolation */
    int r = (int)(rij * t.inv_dx);
    if (r == 0)  ++r;
    const F_FLOAT base = (double)(r+1) * t.dx;
    const F_FLOAT dif = rij - base;

    const cubic_spline_coef vdW = t.d_vdW[r];
    const cubic_spline_coef ele = t.d_ele[r];
    const cubic_spline_coef CEvd = t.d_CEvd[r];
    const cubic_spline_coef CEclmb = t.d_CEclmb[r];

    const F_FLOAT evdwl = ((vdW.d*dif + vdW.c)*dif + vdW.b)*dif +
      vdW.a;

    F_FLOAT ecoul = (((ele.d*dif + ele.c)*dif + ele.b)*dif +
      ele.a)*qi*qj;

    const F_FLOAT fvdwl = ((CEvd.d*dif + CEvd.c)*dif + CEvd.b)*dif +
      CEvd.a;

    F_FLOAT fcoul = (((CEclmb.d*dif+CEclmb.c)*dif+CEclmb.b)*dif +
      CEclmb.a)*qi*qj;

    /* contribution to energy and gradients (atoms and cell)
     * due to geometry-dependent terms in the ACKS2
     * kinetic energy */
    if (acks2_flag) {

      /* kinetic energy terms */
      double xcut = 0.5 * (paramssing(itype).bcut_acks2
                          + paramssing(jtype).bcut_acks2);

      if (rij <= xcut) {
        const F_FLOAT d = rij / xcut;
        const F_FLOAT bond_softness = gp[34] * pow( d, 3.0 )
                                    * pow( 1.0 - d, 6.0 );

        if (bond_softness > 0.0) {
          /* Coulombic energy contribution */
          const F_FLOAT effpot_diff = d_s[NN + i]
                                    - d_s[NN + j];
          const F_FLOAT e_ele = -0.5 * KCALpMOL_to_EV * bond_softness
                                     * SQR( effpot_diff );

          ecoul += e_ele;

          /* forces contribution */
          F_FLOAT d_bond_softness;
          d_bond_softness = gp[34]
                          * 3.0 / xcut * pow( d, 2.0 )
                          * pow( 1.0 - d, 5.0 ) * (1.0 - 3.0 * d);
          d_bond_softness = -0.5 * d_bond_softness
                          * SQR( effpot_diff );
          d_bond_softness = KCALpMOL_to_EV * d_bond_softness
                          / rij;

          fcoul += d_bond_softness;
        }
      }
    }

    const F_FLOAT ftotal = fvdwl + fcoul;
    fxtmp += delx*ftotal;
    fytmp += dely*ftotal;
    fztmp += delz*ftotal;
    a_f(j,0) -= delx*ftotal;
    a_f(j,1) -= dely*ftotal;
    a_f(j,2) -= delz*ftotal;

    if (EVFLAG) {
      if (eflag_global) ev.evdwl += evdwl;
      if (eflag_global) ev.ecoul += ecoul;

      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl+ecoul,-ftotal,delx,dely,delz);
    }
  }

  a_f(i,0) += fxtmp;
  a_f(i,1) += fytmp;
  a_f(i,2) += fztmp;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::allocate_array()
{
  // free scatterview memory
  if (need_dup) {
    dup_dDeltap_self = decltype(dup_dDeltap_self)();
    dup_total_bo     = decltype(dup_total_bo)();
    dup_CdDelta      = decltype(dup_CdDelta)();
  } else {
    ndup_dDeltap_self = decltype(ndup_dDeltap_self)();
    ndup_total_bo     = decltype(ndup_total_bo)();
    ndup_CdDelta      = decltype(ndup_CdDelta)();
  }

  if (cut_hbsq > 0.0) {
    MemKK::realloc_kokkos(d_hb_first,"reaxff/kk:hb_first",nmax);
    MemKK::realloc_kokkos(d_hb_num,"reaxff/kk:hb_num",nmax);
    MemKK::realloc_kokkos(d_hb_list,"reaxff/kk:hb_list",nmax*maxhb);
  }
  MemKK::realloc_kokkos(d_bo_first,"reaxff/kk:bo_first",nmax);
  MemKK::realloc_kokkos(d_bo_num,"reaxff/kk:bo_num",nmax);
  MemKK::realloc_kokkos(d_bo_list,"reaxff/kk:bo_list",nmax*maxbo);

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

  MemKK::realloc_kokkos(d_Cdbo,"reaxff/kk:Cdbo",nmax,3*maxbo);
  MemKK::realloc_kokkos(d_Cdbopi,"reaxff/kk:Cdbopi",nmax,3*maxbo);
  MemKK::realloc_kokkos(d_Cdbopi2,"reaxff/kk:Cdbopi2",nmax,3*maxbo);

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfBlocking<NEIGHFLAG>, const int &ii) const {
  constexpr int blocksize = PairReaxFFKokkos<DeviceType>::build_lists_half_blocksize;

  const auto v_dDeltap_self = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_dDeltap_self),decltype(ndup_dDeltap_self)>::get(dup_dDeltap_self,ndup_dDeltap_self);
  const auto a_dDeltap_self = v_dDeltap_self.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const auto v_total_bo = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_total_bo),decltype(ndup_total_bo)>::get(dup_total_bo,ndup_total_bo);
  const auto a_total_bo = v_total_bo.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3], dBOp_i[3];
  F_FLOAT dDeltap_self_i[3] = {0.0,0.0,0.0};
  F_FLOAT total_bo_i = 0.0;

  int j_index,i_index;
  d_bo_first[i] = i*maxbo;
  const int bo_first_i = d_bo_first[i];

  int ihb = -1;

  int hb_first_i;
  if (cut_hbsq > 0.0) {
    ihb = paramssing(itype).p_hbond;
    if (ihb == 1) {
      d_hb_first[i] = i*maxhb;
      hb_first_i = d_hb_first[i];
    }
  }

  int nnz;
  blocking_t selected_jj[blocksize];
  int jj_current = 0;

  while (jj_current < jnum) {
    nnz = 0;

    while (nnz < blocksize) {
      int jj = jj_current;
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      d_bo_first[j] = j*maxbo;
      d_hb_first[j] = j*maxhb;

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      double cutoffsq;
      if (i < nlocal) cutoffsq = MAX(cut_bosq,cut_hbsq);
      else cutoffsq = cut_bosq;
      if (rsq <= cutoffsq) {
        selected_jj[nnz] = jj_current;
        nnz++;
      }
      jj_current++;

      if (jj_current == jnum) break;
    }

    for (int jj_inner = 0; jj_inner < nnz; jj_inner++) {
      const int jj = selected_jj[jj_inner];
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = type(j);
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // hbond list
      build_hb_list<NEIGHFLAG>(rsq, i, hb_first_i, ihb, j, jtype);

      if (rsq > cut_bosq) continue;

      // bond_list
      const F_FLOAT rij = sqrt(rsq);
      const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
      const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
      const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

      // returns BO_*, C** by reference
      compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
        BO_s, BO_pi, BO_pi2, C12, C34, C56);

      BO = BO_s + BO_pi + BO_pi2;
      if (BO < bo_cut) continue;

      int ii_index = -1;
      int jj_index = -1;
      if (build_bo_list<NEIGHFLAG>(bo_first_i, i, j, i_index, j_index, ii_index, jj_index)) {

        // from BondOrder1

        d_BO(i,jj_index) = BO;
        d_BO_s(i,jj_index) = BO_s;

        d_BO(j,ii_index) = BO;
        d_BO_s(j,ii_index) = BO_s;

        d_BO_pi(j,ii_index) = BO_pi;
        d_BO_pi2(j,ii_index) = BO_pi2;

        d_BO_pi(i,jj_index) = BO_pi;
        d_BO_pi2(i,jj_index) = BO_pi2;

        F_FLOAT Cln_BOp_s = p_bo2 * C12 / rij / rij;
        F_FLOAT Cln_BOp_pi = p_bo4 * C34 / rij / rij;
        F_FLOAT Cln_BOp_pi2 = p_bo6 * C56 / rij / rij;

        if (nlocal == 0)
          Cln_BOp_s = Cln_BOp_pi = Cln_BOp_pi2 = 0.0;

        for (int d = 0; d < 3; d++) dBOp_i[d] = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2)*delij[d];
        for (int d = 0; d < 3; d++) dDeltap_self_i[d] += dBOp_i[d];
        for (int d = 0; d < 3; d++) a_dDeltap_self(j,d) += -dBOp_i[d];

        d_dln_BOp_pi(i,jj_index) = -(BO_pi*Cln_BOp_pi);
        d_dln_BOp_pi(j,ii_index) = -(BO_pi*Cln_BOp_pi);

        d_dln_BOp_pi2(i,jj_index) = -(BO_pi2*Cln_BOp_pi2);
        d_dln_BOp_pi2(j,ii_index) = -(BO_pi2*Cln_BOp_pi2);

        d_dBOp(i,jj_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);
        d_dBOp(j,ii_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);
        d_BO(i,jj_index) = BO - bo_cut;
        d_BO(j,ii_index) = BO - bo_cut;
        d_BO_s(i,jj_index) = BO_s - bo_cut;
        d_BO_s(j,ii_index) = BO_s - bo_cut;
        total_bo_i += (BO - bo_cut);
        a_total_bo[j] += (BO - bo_cut);
      }
    }
  }

  for (int d = 0; d < 3; d++)
    a_dDeltap_self(i,d) += dDeltap_self_i[d];

  a_total_bo[i] += total_bo_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfBlockingPreview<NEIGHFLAG>, const int &ii) const {
  constexpr int blocksize = PairReaxFFKokkos<DeviceType>::build_lists_half_blocksize;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3];

  int j_index,i_index;
  d_bo_first[i] = i*maxbo;
  const int bo_first_i = d_bo_first[i];

  int ihb = -1;

  int hb_first_i;
  if (cut_hbsq > 0.0) {
    ihb = paramssing(itype).p_hbond;
    if (ihb == 1) {
      d_hb_first[i] = i*maxhb;
      hb_first_i = d_hb_first[i];
    }
  }

  int nnz;
  blocking_t selected_jj[blocksize];
  int jj_current = 0;

  double cutoffsq;
  if (i < nlocal) cutoffsq = MAX(cut_bosq,cut_hbsq);
  else cutoffsq = cut_bosq;

  while (jj_current < jnum) {
    nnz = 0;

    while (nnz < blocksize) {
      int jj = jj_current;
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      d_bo_first[j] = j*maxbo;
      d_hb_first[j] = j*maxhb;

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      if (rsq <= cutoffsq) {
        selected_jj[nnz] = jj_current;
        nnz++;
      }
      jj_current++;

      if (jj_current == jnum) break;
    }

    for (int jj_inner = 0; jj_inner < nnz; jj_inner++) {
      const int jj = selected_jj[jj_inner];
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = type(j);
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // hbond list
      build_hb_list<NEIGHFLAG>(rsq, i, hb_first_i, ihb, j, jtype);

      if (rsq > cut_bosq) continue;

      // bond_list
      const F_FLOAT rij = sqrt(rsq);
      const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
      const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
      const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

      // returns BO_*, C** by reference
      compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
        BO_s, BO_pi, BO_pi2, C12, C34, C56);

      BO = BO_s + BO_pi + BO_pi2;
      if (BO < bo_cut) continue;

      int ii_index = -1;
      int jj_index = -1;
      build_bo_list<NEIGHFLAG>(bo_first_i, i, j, i_index, j_index, ii_index, jj_index);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfPreview<NEIGHFLAG>, const int &ii) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3];

  int j_index,i_index;
  d_bo_first[i] = i*maxbo;
  const int bo_first_i = d_bo_first[i];

  int ihb = -1;

  int hb_first_i;
  if (cut_hbsq > 0.0) {
    ihb = paramssing(itype).p_hbond;
    if (ihb == 1) {
      d_hb_first[i] = i*maxhb;
      hb_first_i = d_hb_first[i];
    }
  }

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    d_bo_first[j] = j*maxbo;
    d_hb_first[j] = j*maxhb;
    const int jtype = type(j);

    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;
    const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

    // hbond list
    build_hb_list<NEIGHFLAG>(rsq, i, hb_first_i, ihb, j, jtype);

    if (rsq > cut_bosq) continue;

    // bond_list
    const F_FLOAT rij = sqrt(rsq);
    const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
    const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
    const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

    // returns BO_*, C** by reference
    compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
      BO_s, BO_pi, BO_pi2, C12, C34, C56);

    BO = BO_s + BO_pi + BO_pi2;
    if (BO < bo_cut) continue;

    int ii_index = -1;
    int jj_index = -1;

    build_bo_list<NEIGHFLAG>(bo_first_i, i, j, i_index, j_index, ii_index, jj_index);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::build_hb_list(F_FLOAT rsq, int i, int hb_first_i, int ihb, int j, int jtype) const {

  int i_index, j_index;
  int jhb = -1;
  if (i < nlocal && cut_hbsq > 0.0 && (ihb == 1 || ihb == 2) && rsq <= cut_hbsq) {
    jhb = paramssing(jtype).p_hbond;
    if (ihb == 1 && jhb == 2) {
      if (NEIGHFLAG == HALF) {
        j_index = hb_first_i + d_hb_num[i];
        d_hb_num[i]++;
      } else
        j_index = hb_first_i + Kokkos::atomic_fetch_add(&d_hb_num[i],1);

      const int jj_index = j_index - hb_first_i;

      if (jj_index >= maxhb)
        d_resize_hb() = MAX(d_resize_hb(),jj_index+1);
      else
        d_hb_list[j_index] = j;
    } else if (j < nlocal && ihb == 2 && jhb == 1) {
      if (NEIGHFLAG == HALF) {
        i_index = d_hb_first[j] + d_hb_num[j];
        d_hb_num[j]++;
      } else
        i_index = d_hb_first[j] + Kokkos::atomic_fetch_add(&d_hb_num[j],1);

      const int ii_index = i_index - d_hb_first[j];

      if (ii_index >= maxhb)
        d_resize_hb() = MAX(d_resize_hb(),ii_index+1);
      else
        d_hb_list[i_index] = i;
    }
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
bool PairReaxFFKokkos<DeviceType>::build_bo_list(int bo_first_i, int i, int j, int i_index, int j_index, int& ii_index, int& jj_index) const {

  if (NEIGHFLAG == HALF) {
    j_index = bo_first_i + d_bo_num[i];
    i_index = d_bo_first[j] + d_bo_num[j];
    d_bo_num[i]++;
    d_bo_num[j]++;
  } else {
    j_index = bo_first_i + Kokkos::atomic_fetch_add(&d_bo_num[i],1);
    i_index = d_bo_first[j] + Kokkos::atomic_fetch_add(&d_bo_num[j],1);
  }

  jj_index = j_index - bo_first_i;
  ii_index = i_index - d_bo_first[j];

  bool set_dB_flag = true;

  if (jj_index >= maxbo || ii_index >= maxbo) {
    const int max_val = MAX(ii_index+1,jj_index+1);
    d_resize_bo() = MAX(d_resize_bo(),max_val);
    set_dB_flag = false;
  } else {
    d_bo_list[j_index] = j;
    d_bo_list[i_index] = i;
    set_dB_flag = true;
  }

  return set_dB_flag;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsFull, const int &ii) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3], dBOp_i[3];
  F_FLOAT dDeltap_self_i[3] = {0.0,0.0,0.0};
  F_FLOAT total_bo_i = 0.0;

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];
  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int jtype = type(j);
    const int j_index = jj - j_start;
    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;
    const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
    const F_FLOAT rsq_inv = 1.0 / rsq;

    // bond_list
    const F_FLOAT rij = sqrt(rsq);
    const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
    const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
    const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

    // returns BO_*, C** by reference
    compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
      BO_s, BO_pi, BO_pi2, C12, C34, C56);

    BO = BO_s + BO_pi + BO_pi2;

    // from BondOrder1

    d_BO(i,j_index) = BO;
    d_BO_s(i,j_index) = BO_s;
    d_BO_pi(i,j_index) = BO_pi;
    d_BO_pi2(i,j_index) = BO_pi2;

    F_FLOAT Cln_BOp_s = p_bo2 * C12 * rsq_inv;
    F_FLOAT Cln_BOp_pi = p_bo4 * C34 * rsq_inv;
    F_FLOAT Cln_BOp_pi2 = p_bo6 * C56 * rsq_inv;

    if (nlocal == 0)
      Cln_BOp_s = Cln_BOp_pi = Cln_BOp_pi2 = 0.0;

    for (int d = 0; d < 3; d++) dBOp_i[d] = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2)*delij[d];
    for (int d = 0; d < 3; d++) dDeltap_self_i[d] += dBOp_i[d];


    d_dln_BOp_pi(i,j_index) = -(BO_pi*Cln_BOp_pi);
    d_dln_BOp_pi2(i,j_index) = -(BO_pi2*Cln_BOp_pi2);
    d_dBOp(i,j_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);

    d_BO(i,j_index) = BO - bo_cut;
    d_BO_s(i,j_index) = BO_s - bo_cut;
    total_bo_i += (BO - bo_cut);
  }

  for (int d = 0; d < 3; d++)
    d_dDeltap_self(i,d) = dDeltap_self_i[d];

  d_total_bo[i] = total_bo_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::compute_bo(F_FLOAT rij, int itype, int jtype, F_FLOAT p_bo2, F_FLOAT p_bo4, F_FLOAT p_bo6,
  F_FLOAT& BO_s, F_FLOAT& BO_pi, F_FLOAT& BO_pi2, F_FLOAT& C12, F_FLOAT& C34, F_FLOAT& C56) const {

  const F_FLOAT p_bo1 = paramstwbp(itype,jtype).p_bo1;
  const F_FLOAT p_bo3 = paramstwbp(itype,jtype).p_bo3;
  const F_FLOAT p_bo5 = paramstwbp(itype,jtype).p_bo5;
  const F_FLOAT r_s = paramstwbp(itype,jtype).r_s;
  const F_FLOAT r_pi = paramstwbp(itype,jtype).r_pi;
  const F_FLOAT r_pi2 = paramstwbp(itype,jtype).r_pi2;

  if (paramssing(itype).r_s > 0.0  && paramssing(jtype).r_s > 0.0) {
    C12 = p_bo1 * ((p_bo2 != 0) ? (pow(rij/r_s,p_bo2)) : 1.0);
    BO_s = (1.0+bo_cut)*exp(C12);
  } else BO_s = C12 = 0.0;

  if (paramssing(itype).r_pi > 0.0  && paramssing(jtype).r_pi > 0.0) {
    C34 = p_bo3 * ((p_bo4 != 0) ? (pow(rij/r_pi,p_bo4)) : 1.0);
    BO_pi = exp(C34);
  } else BO_pi = C34 = 0.0;

  if (paramssing(itype).r_pi2 > 0.0  && paramssing(jtype).r_pi2 > 0.0) {
    C56 = p_bo5 * ((p_bo6 != 0) ? (pow(rij/r_pi2,p_bo6)) : 1.0);
    BO_pi2 = exp(C56);
  } else BO_pi2 = C56 = 0.0;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder1, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);

  const F_FLOAT val_i = paramssing(itype).valency;
  d_Deltap[i] = d_total_bo[i] - val_i;
  d_Deltap_boc[i] = d_total_bo[i] - paramssing(itype).valency_val;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder2, const int &ii) const {

  F_FLOAT exp_p1i, exp_p2i, exp_p1j, exp_p2j, f1, f2, f3, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  F_FLOAT f4, f5, exp_f4, exp_f5, f4f5, Cf45_ij, Cf45_ji;
  F_FLOAT A0_ij, A1_ij, A2_ij, A3_ij, A2_ji, A3_ji;

  const int i = d_ilist[ii];
  const int itype = type(i);
  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  const F_FLOAT val_i = paramssing(itype).valency;

  d_total_bo[i] = 0.0;
  F_FLOAT total_bo = 0.0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int jtype = type(j);
    const int j_index = jj - j_start;
    const int i_index = maxbo+j_index;

    // calculate corrected BO and total bond order

    const F_FLOAT val_j = paramssing(jtype).valency;
    const F_FLOAT ovc = paramstwbp(itype,jtype).ovc;
    const F_FLOAT v13cor = paramstwbp(itype,jtype).v13cor;
    const F_FLOAT p_boc3 = paramstwbp(itype,jtype).p_boc3;
    const F_FLOAT p_boc4 = paramstwbp(itype,jtype).p_boc4;
    const F_FLOAT p_boc5 = paramstwbp(itype,jtype).p_boc5;

    if (ovc < 0.001 && v13cor < 0.001) {
      d_C1dbo(i,j_index) = 1.0;
      d_C2dbo(i,j_index) = 0.0;
      d_C3dbo(i,j_index) = 0.0;
      d_C1dbopi(i,j_index) = 1.0;
      d_C2dbopi(i,j_index) = 0.0;
      d_C3dbopi(i,j_index) = 0.0;
      d_C4dbopi(i,j_index) = 0.0;
      d_C1dbopi2(i,j_index) = 1.0;
      d_C2dbopi2(i,j_index) = 0.0;
      d_C3dbopi2(i,j_index) = 0.0;
      d_C4dbopi2(i,j_index) = 0.0;
    } else {
      if (ovc >= 0.001) {
        exp_p1i = exp(-p_boc1 * d_Deltap[i]);
        exp_p2i = exp(-p_boc2 * d_Deltap[i]);
        exp_p1j = exp(-p_boc1 * d_Deltap[j]);
        exp_p2j = exp(-p_boc2 * d_Deltap[j]);

        f2 = exp_p1i + exp_p1j;
        f3 = -1.0/p_boc2*log(0.5*(exp_p2i+exp_p2j));
        f1 = 0.5 * ((val_i + f2)/(val_i + f2 + f3) + (val_j + f2)/(val_j + f2 + f3));
        u1_ij = val_i + f2 + f3;
        u1_ji = val_j + f2 + f3;
        Cf1A_ij = 0.5 * f3 * (1.0/(u1_ij*u1_ij)+1.0/(u1_ji*u1_ji));
        Cf1B_ij = -0.5 * ((u1_ij - f3)/(u1_ij*u1_ij)+(u1_ji - f3)/(u1_ji*u1_ji));
        Cf1_ij = 0.5 * (-p_boc1 * exp_p1i / u1_ij - ((val_i+f2) / (u1_ij*u1_ij)) *
                       (-p_boc1 * exp_p1i + exp_p2i / (exp_p2i + exp_p2j)) +
                        -p_boc1 * exp_p1i / u1_ji - ((val_j+f2) / (u1_ji*u1_ji)) *
                       (-p_boc1 * exp_p1i + exp_p2i / (exp_p2i + exp_p2j)));
        Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j + Cf1B_ij * exp_p2j / (exp_p2i + exp_p2j);
      } else {
        f1 = 1.0;
        Cf1_ij = Cf1_ji = 0.0;
      }

      if (v13cor >= 0.001) {
        exp_f4 =exp(-(p_boc4*(d_BO(i,j_index)*d_BO(i,j_index))-d_Deltap_boc[i])*p_boc3+p_boc5);
        exp_f5 =exp(-(p_boc4*(d_BO(i,j_index)*d_BO(i,j_index))-d_Deltap_boc[j])*p_boc3+p_boc5);
        f4 = 1. / (1. + exp_f4);
        f5 = 1. / (1. + exp_f5);
        f4f5 = f4 * f5;

        Cf45_ij = -f4 * exp_f4;
        Cf45_ji = -f5 * exp_f5;
      } else {
        f4 = f5 = f4f5 = 1.0;
        Cf45_ij = Cf45_ji = 0.0;
      }

      A0_ij = f1 * f4f5;
      A1_ij = -2 * p_boc3 * p_boc4 * d_BO(i,j_index) * (Cf45_ij + Cf45_ji);
      A2_ij = Cf1_ij / f1 + p_boc3 * Cf45_ij;
      A2_ji = Cf1_ji / f1 + p_boc3 * Cf45_ji;
      A3_ij = A2_ij + Cf1_ij / f1;
      A3_ji = A2_ji + Cf1_ji / f1;

      d_BO(i,j_index) = d_BO(i,j_index) * A0_ij;
      d_BO_pi(i,j_index) = d_BO_pi(i,j_index) * A0_ij * f1;
      d_BO_pi2(i,j_index) = d_BO_pi2(i,j_index) * A0_ij * f1;
      d_BO_s(i,j_index) = d_BO(i,j_index)-(d_BO_pi(i,j_index)+d_BO_pi2(i,j_index));

      d_C1dbo(i,j_index) = A0_ij + d_BO(i,j_index) * A1_ij;
      d_C2dbo(i,j_index) = d_BO(i,j_index) * A2_ij;
      d_C3dbo(i,j_index) = d_BO(i,j_index) * A2_ji;

      d_C1dbopi(i,j_index) = f1*f1*f4*f5;
      d_C2dbopi(i,j_index) = d_BO_pi(i,j_index) * A1_ij;
      d_C3dbopi(i,j_index) = d_BO_pi(i,j_index) * A3_ij;
      d_C4dbopi(i,j_index) = d_BO_pi(i,j_index) * A3_ji;

      d_C1dbopi2(i,j_index) = f1*f1*f4*f5;
      d_C2dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A1_ij;
      d_C3dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A3_ij;
      d_C4dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A3_ji;
    }

    if (d_BO(i,j_index) < 1e-10) d_BO(i,j_index) = 0.0;
    if (d_BO_s(i,j_index) < 1e-10) d_BO_s(i,j_index) = 0.0;
    if (d_BO_pi(i,j_index) < 1e-10) d_BO_pi(i,j_index) = 0.0;
    if (d_BO_pi2(i,j_index) < 1e-10) d_BO_pi2(i,j_index) = 0.0;

    total_bo += d_BO(i,j_index);

    d_Cdbo(i,j_index) = 0.0;
    d_Cdbopi(i,j_index) = 0.0;
    d_Cdbopi2(i,j_index) = 0.0;
    d_Cdbo(j,i_index) = 0.0;
    d_Cdbopi(j,i_index) = 0.0;
    d_Cdbopi2(j,i_index) = 0.0;

    d_CdDelta[j] = 0.0;
  }
  d_CdDelta[i] = 0.0;
  d_total_bo[i] += total_bo;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder3, const int &ii) const {
  // bot part of BO()

  const int i = d_ilist[ii];
  const int itype = type(i);
  F_FLOAT nlp_temp;

  d_Delta[i] = d_total_bo[i] - paramssing(itype).valency;
  const F_FLOAT Delta_e = d_total_bo[i] - paramssing(itype).valency_e;
  d_Delta_boc[i] = d_total_bo[i] - paramssing(itype).valency_boc;

  const F_FLOAT vlpex = Delta_e - 2.0 * (int)(Delta_e/2.0);
  const F_FLOAT explp1 = exp(-gp[15] * SQR(2.0 + vlpex));
  const F_FLOAT nlp = explp1 - (int)(Delta_e / 2.0);
  d_Delta_lp[i] = paramssing(itype).nlp_opt - nlp;
  const F_FLOAT Clp = 2.0 * gp[15] * explp1 * (2.0 + vlpex);
  d_dDelta_lp[i] = Clp;

  if (paramssing(itype).mass > 21.0) {
    nlp_temp = 0.5 * (paramssing(itype).valency_e - paramssing(itype).valency);
    d_Delta_lp_temp[i] = paramssing(itype).nlp_opt - nlp_temp;
  } else {
    nlp_temp = nlp;
    d_Delta_lp_temp[i] = paramssing(itype).nlp_opt - nlp_temp;
  }

  d_sum_ovun(i,1) = 0.0;
  d_sum_ovun(i,2) = 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti1, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT imass = paramssing(itype).mass;
  F_FLOAT dfvl;

  if (imass > 21.0) dfvl = 0.0;
  else dfvl = 1.0;

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  F_FLOAT sum_ovun1 = 0.0;
  F_FLOAT sum_ovun2 = 0.0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int jtype = type(j);
    const int j_index = jj - j_start;

    sum_ovun1 += paramstwbp(itype,jtype).p_ovun1 * paramstwbp(itype,jtype).De_s * d_BO(i,j_index);
    sum_ovun2 += (d_Delta[j] - dfvl * d_Delta_lp_temp[j]) * (d_BO_pi(i,j_index) + d_BO_pi2(i,j_index));
  }
  d_sum_ovun(i,1) += sum_ovun1;
  d_sum_ovun(i,2) += sum_ovun2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT imass = paramssing(itype).mass;
  const F_FLOAT val_i = paramssing(itype).valency;

  F_FLOAT dfvl;
  if (imass > 21.0) dfvl = 0.0;
  else dfvl = 1.0;

  F_FLOAT e_lp, e_ov, e_un;
  F_FLOAT CEover1, CEover2, CEover3, CEover4;
  F_FLOAT CEunder1, CEunder2, CEunder3, CEunder4;
  const F_FLOAT p_lp3 = gp[5];
  const F_FLOAT p_ovun2 = paramssing(itype).p_ovun2;
  const F_FLOAT p_ovun3 = gp[32];
  const F_FLOAT p_ovun4 = gp[31];
  const F_FLOAT p_ovun5 = paramssing(itype).p_ovun5;
  const F_FLOAT p_ovun6 = gp[6];
  const F_FLOAT p_ovun7 = gp[8];
  const F_FLOAT p_ovun8 = gp[9];

  // lone pair
  const F_FLOAT p_lp2 = paramssing(itype).p_lp2;
  const F_FLOAT expvd2 = exp(-75 * d_Delta_lp[i]);
  const F_FLOAT inv_expvd2 = 1.0 / (1.0+expvd2);

  int numbonds = d_bo_num[i];

  e_lp = 0.0;
  if (numbonds > 0 || enobondsflag)
    e_lp = p_lp2 * d_Delta_lp[i] * inv_expvd2;
  const F_FLOAT dElp = p_lp2 * inv_expvd2 + 75.0 * p_lp2 * d_Delta_lp[i] * expvd2 * inv_expvd2*inv_expvd2;
  const F_FLOAT CElp = dElp * d_dDelta_lp[i];

  if (numbonds > 0 || enobondsflag)
    a_CdDelta[i] += CElp;

  if (EFLAG && eflag_global) ev.ereax[0] += e_lp;
  //if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_lp,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_lp);

  // over coordination
  const F_FLOAT exp_ovun1 = p_ovun3 * exp(p_ovun4 * d_sum_ovun(i,2));
  const F_FLOAT inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
  const F_FLOAT Delta_lpcorr  = d_Delta[i] - (dfvl * d_Delta_lp_temp[i]) * inv_exp_ovun1;

  const F_FLOAT exp_ovun2 = exp(p_ovun2 * Delta_lpcorr);
  const F_FLOAT inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);
  const F_FLOAT DlpVi = 1.0 / (Delta_lpcorr + val_i + 1e-8);

  CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;
  e_ov = d_sum_ovun(i,1) * CEover1;

  if (EFLAG && eflag_global) ev.ereax[1] += e_ov;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_ov,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_ov);

  CEover2 = d_sum_ovun(i,1) * DlpVi * inv_exp_ovun2 *
    (1.0 - Delta_lpcorr * (DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2));
  CEover3 = CEover2 * (1.0 - dfvl * d_dDelta_lp[i] * inv_exp_ovun1);
  CEover4 = CEover2 * (dfvl * d_Delta_lp_temp[i]) * p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);

  // under coordination

  const F_FLOAT exp_ovun2n = 1.0 / exp_ovun2;
  const F_FLOAT exp_ovun6 = exp(p_ovun6 * Delta_lpcorr);
  const F_FLOAT exp_ovun8 = p_ovun7 * exp(p_ovun8 * d_sum_ovun(i,2));
  const F_FLOAT inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
  const F_FLOAT inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

  e_un = 0;
  if (numbonds > 0 || enobondsflag)
    e_un = -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;

  if (EFLAG && eflag_global) ev.ereax[2] += e_un;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_un,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_un);

  CEunder1 = inv_exp_ovun2n *
    (p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 + p_ovun2 * e_un * exp_ovun2n);
  CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
  CEunder3 = CEunder1 * (1.0 - dfvl * d_dDelta_lp[i] * inv_exp_ovun1);
  CEunder4 = CEunder1 * (dfvl * d_Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * inv_exp_ovun1 * inv_exp_ovun1 + CEunder2;

  const F_FLOAT eng_tmp = e_lp + e_ov + e_un;
  if (eflag_atom) this->template e_tally_single<NEIGHFLAG>(ev,i,eng_tmp);

  // multibody forces

  a_CdDelta[i] += CEover3;
  if (numbonds > 0 || enobondsflag)
    a_CdDelta[i] += CEunder3;

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  F_FLOAT CdDelta_i = 0.0;
  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int jtype = type(j);
    const F_FLOAT jmass = paramssing(jtype).mass;
    const int j_index = jj - j_start;
    const F_FLOAT De_s = paramstwbp(itype,jtype).De_s;

    // multibody lone pair: correction for C2
    if (p_lp3 > 0.001 && imass == 12.0 && jmass == 12.0) {
      const F_FLOAT Di = d_Delta[i];
      const F_FLOAT vov3 = d_BO(i,j_index) - Di - 0.040*pow(Di,4.0);
      if (vov3 > 3.0) {
        const F_FLOAT e_lph = p_lp3 * (vov3-3.0)*(vov3-3.0);
        const F_FLOAT deahu2dbo = 2.0 * p_lp3 * (vov3 - 3.0);
        const F_FLOAT deahu2dsbo = 2.0 * p_lp3 * (vov3 - 3.0) * (-1.0 - 0.16 * pow(Di,3.0));
        d_Cdbo(i,j_index) += deahu2dbo;
        CdDelta_i += deahu2dsbo;

        if (EFLAG) {
          if (eflag_global) ev.ereax[0] += e_lph;
          if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,e_lph);
        }
      }
    }

    // over/under coordination forces merged together
    const F_FLOAT p_ovun1 = paramstwbp(itype,jtype).p_ovun1;
    a_CdDelta[j] += (CEover4 + CEunder4) * (1.0 - dfvl * d_dDelta_lp[j]) * (d_BO_pi(i,j_index) + d_BO_pi2(i,j_index));
    d_Cdbo(i,j_index) += CEover1 * p_ovun1 * De_s;
    d_Cdbopi(i,j_index) += (CEover4 + CEunder4) * (d_Delta[j] - dfvl*d_Delta_lp_temp[j]);
    d_Cdbopi2(i,j_index) += (CEover4 + CEunder4) * (d_Delta[j] - dfvl*d_Delta_lp_temp[j]);
  }
  a_CdDelta[i] += CdDelta_i;

}

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EFLAG>(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxCountAngularTorsion<POPULATE>, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  if (POPULATE) {
    // Computes and stores SBO2, CSBO2, dSBO1, dSBO2
    compute_angular_sbo(i, itype, j_start, j_end);
  }

  // Angular

  // Count buffer size for `i`
  int location_angular = 0; // dummy declaration
  int count_angular = preprocess_angular<false>(i, itype, j_start, j_end, location_angular);
  location_angular = Kokkos::atomic_fetch_add(&d_count_angular_torsion(0), count_angular);

  if (POPULATE) {
    // Fill buffer for `i`
    preprocess_angular<true>(i, itype, j_start, j_end, location_angular);
  }

  // Torsion

  const tagint itag = tag(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  // Count buffer size for `i`
  int location_torsion = 0; // dummy declaration
  int count_torsion = preprocess_torsion<false>(i, itype, itag, xtmp, ytmp, ztmp, j_start, j_end, location_torsion);
  location_torsion = Kokkos::atomic_fetch_add(&d_count_angular_torsion(1), count_torsion);

  if (POPULATE) {
    // Fill buffer for `i`
    preprocess_torsion<true>(i, itype, itag, xtmp, ytmp, ztmp, j_start, j_end, location_torsion);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::compute_angular_sbo(int i, int itype, int j_start, int j_end) const {

  F_FLOAT SBO2, CSBO2, dSBO1, dSBO2;

  const F_FLOAT p_val8 = gp[33];
  const F_FLOAT p_val9 = gp[16];

  F_FLOAT SBOp = 0.0;
  F_FLOAT prod_SBO = 1.0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int j_index = jj - j_start;
    const F_FLOAT bo_ij = d_BO(i,j_index);

    SBOp += (d_BO_pi(i,j_index) + d_BO_pi2(i,j_index));
    F_FLOAT temp = SQR(bo_ij);
    temp *= temp;
    temp *= temp;
    prod_SBO *= exp(-temp);
  }

  F_FLOAT vlpadj;

  const F_FLOAT Delta_e = d_total_bo[i] - paramssing(itype).valency_e;
  const F_FLOAT vlpex = Delta_e - 2.0 * (int)(Delta_e/2.0);
  const F_FLOAT explp1 = exp(-gp[15] * SQR(2.0 + vlpex));
  const F_FLOAT nlp = explp1 - (int)(Delta_e / 2.0);
  if (vlpex >= 0.0) {
    vlpadj = 0.0;
    dSBO2 = prod_SBO - 1.0;
  } else {
    vlpadj = nlp;
    dSBO2 = (prod_SBO - 1.0) * (1.0 - p_val8 * d_dDelta_lp[i]);
  }

  const F_FLOAT SBO = SBOp + (1.0 - prod_SBO) * (-d_Delta_boc[i] - p_val8 * vlpadj);
  dSBO1 = -8.0 * prod_SBO * (d_Delta_boc[i] + p_val8 * vlpadj);

  if (SBO <= 0.0) {
    SBO2 = 0.0;
    CSBO2 = 0.0;
  } else if (SBO > 0.0 && SBO <= 1.0) {
    SBO2 = pow(SBO, p_val9);
    CSBO2 = p_val9 * pow(SBO, p_val9 - 1.0);
  } else if (SBO > 1.0 && SBO < 2.0) {
    SBO2 = 2.0 - pow(2.0-SBO, p_val9);
    CSBO2 = p_val9 * pow(2.0 - SBO, p_val9 - 1.0);
  } else {
    SBO2 = 2.0;
    CSBO2 = 0.0;
  }

  d_angular_intermediates(i,0) = SBO2;
  d_angular_intermediates(i,1) = CSBO2;
  d_angular_intermediates(i,2) = dSBO1;
  d_angular_intermediates(i,3) = dSBO2;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
int PairReaxFFKokkos<DeviceType>::preprocess_angular(int i, int itype, int j_start, int j_end, int location_angular) const {

  int count_angular = 0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int j_index = jj - j_start;
    const F_FLOAT bo_ij = d_BO(i,j_index);

    if (bo_ij <= thb_cut) continue;
    if (i >= nlocal && j >= nlocal) continue;

    const int i_index = maxbo + j_index;
    const int jtype = type(j);

    for (int kk = jj+1; kk < j_end; kk++) {
    //for (int kk = j_start; kk < j_end; kk++) {
      int k = d_bo_list[kk];
      k &= NEIGHMASK;
      if (k == j) continue;

      const int k_index = kk - j_start;
      const F_FLOAT bo_ik = d_BO(i,k_index);

      if (bo_ij <= thb_cut || bo_ik <= thb_cut || bo_ij * bo_ik <= thb_cutsq) continue;

      const int ktype = type(k);

      F_FLOAT p_val1 = paramsthbp(jtype,itype,ktype).p_val1;

      if (fabs(p_val1) <= 0.001) continue;

      if (POPULATE) {
        reax_int4 pack;

        // First pack stores i, j, k, and j_start
        pack.i0 = i;
        pack.i1 = j;
        pack.i2 = k;
        pack.i3 = j_start;
        d_angular_pack(location_angular, 0) = pack;

        // Second pack stores i_index, j_index, k_index, and j_end
        pack.i0 = i_index;
        pack.i1 = j_index;
        pack.i2 = k_index;
        pack.i3 = j_end;
        d_angular_pack(location_angular, 1) = pack;

        location_angular++;
      } else {
        count_angular++;
      }
    }
  }

  return count_angular;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
int PairReaxFFKokkos<DeviceType>::preprocess_torsion(int i, int /*itype*/, tagint itag,
  F_FLOAT xtmp, F_FLOAT ytmp, F_FLOAT ztmp, int j_start, int j_end, int location_torsion) const {

  // in reaxff_torsion_angles: j = i, k = j, i = k;

  int count_torsion = 0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const tagint jtag = tag(j);
    const int j_index = jj - j_start;

    // skip half of the interactions
    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const F_FLOAT bo_ij = d_BO(i,j_index);
    if (bo_ij < thb_cut) continue;

    const int l_start = d_bo_first[j];
    const int l_end = l_start + d_bo_num[j];

    for (int kk = j_start; kk < j_end; kk++) {
      int k = d_bo_list[kk];
      k &= NEIGHMASK;
      if (k == j) continue;
      const int k_index = kk - j_start;

      const F_FLOAT bo_ik = d_BO(i,k_index);
      if (bo_ik < thb_cut) continue;

      for (int ll = l_start; ll < l_end; ll++) {
        int l = d_bo_list[ll];
        l &= NEIGHMASK;
        if (l == i) continue;
        const int l_index = ll - l_start;

        const F_FLOAT bo_jl = d_BO(j,l_index);
        if (l == k || bo_jl < thb_cut || bo_ij*bo_ik*bo_jl < thb_cut) continue;

        if (POPULATE) {
          reax_int4 pack;

          pack.i0 = i;
          pack.i1 = j;
          pack.i2 = k;
          pack.i3 = l;
          d_torsion_pack(location_torsion, 0) = pack;

          pack.i0 = 0; // no i_index
          pack.i1 = j_index;
          pack.i2 = k_index;
          pack.i3 = l_index;
          d_torsion_pack(location_torsion, 1) = pack;

          location_torsion++;
        } else {
          count_torsion++;
        }
      }
    }
  }

  return count_torsion;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int &apack, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi2 = d_Cdbopi2;

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT temp, temp_bo_jt, pBOjt7;
  F_FLOAT p_val1, p_val2, p_val3, p_val4, p_val5;
  F_FLOAT p_val6, p_val7, p_val10;
  F_FLOAT p_pen1, p_pen2, p_pen3, p_pen4;
  F_FLOAT p_coa1, p_coa2, p_coa3, p_coa4;
  F_FLOAT trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
  F_FLOAT exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
  F_FLOAT dSBO1, dSBO2, SBO2, CSBO2;
  F_FLOAT CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
  F_FLOAT CEpen1, CEpen2, CEpen3;
  F_FLOAT e_ang, e_coa, e_pen;
  F_FLOAT CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
  F_FLOAT Cf7ij, Cf7jk, Cf8j, Cf9j;
  F_FLOAT f7_ij, f7_jk, f8_Dj, f9_Dj;
  F_FLOAT Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta;
  F_FLOAT BOA_ij, BOA_ik, rij, bo_ij, bo_ik;
  F_FLOAT dcos_theta_di[3], dcos_theta_dj[3], dcos_theta_dk[3];
  F_FLOAT eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
  F_FLOAT delij[3], delik[3], delji[3], delki[3];

  p_val6 = gp[14];
  p_val10 = gp[17];

  p_pen2 = gp[19];
  p_pen3 = gp[20];
  p_pen4 = gp[21];

  p_coa2 = gp[2];
  p_coa3 = gp[38];
  p_coa4 = gp[30];

  reax_int4 pack = d_angular_pack(apack,0);
  const int i = pack.i0;
  const int j = pack.i1;
  const int k = pack.i2;
  const int j_start = pack.i3;

  pack = d_angular_pack(apack, 1);
  const int i_index = pack.i0;
  const int j_index = pack.i1;
  const int k_index = pack.i2;
  const int j_end = pack.i3;

  const int itype = type(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  p_val3 = paramssing(itype).p_val3;
  p_val5 = paramssing(itype).p_val5;

  const F_FLOAT Delta_val = d_total_bo[i] - paramssing(itype).valency_val;

  SBO2 = d_angular_intermediates(i, 0);
  CSBO2 = d_angular_intermediates(i, 1);
  dSBO1 = d_angular_intermediates(i, 2);
  dSBO2 = d_angular_intermediates(i, 3);

  expval6 = exp(p_val6 * d_Delta_boc[i]);

  F_FLOAT CdDelta_i = 0.0;
  F_FLOAT fitmp[3],fjtmp[3];
  for (int j = 0; j < 3; j++) fitmp[j] = 0.0;

  delij[0] = x(j,0) - xtmp;
  delij[1] = x(j,1) - ytmp;
  delij[2] = x(j,2) - ztmp;
  const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
  rij = sqrt(rsqij);
  bo_ij = d_BO(i,j_index);

  BOA_ij = bo_ij - thb_cut;

  const int jtype = type(j);

  F_FLOAT CdDelta_j = 0.0;
  for (int k = 0; k < 3; k++) fjtmp[k] = 0.0;

  delik[0] = x(k,0) - xtmp;
  delik[1] = x(k,1) - ytmp;
  delik[2] = x(k,2) - ztmp;
  const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
  const F_FLOAT rik = sqrt(rsqik);
  bo_ik = d_BO(i,k_index);
  BOA_ik   = bo_ik - thb_cut;

  const int ktype = type(k);

  // theta and derivatives

  cos_theta = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
  if (cos_theta > 1.0) cos_theta  = 1.0;
  if (cos_theta < -1.0) cos_theta  = -1.0;
  theta = acos(cos_theta);

  const F_FLOAT inv_dists = 1.0 / (rij * rik);
  const F_FLOAT Cdot_inv3 = cos_theta * inv_dists * inv_dists;

  for (int t = 0; t < 3; t++) {
    dcos_theta_di[t] = -(delik[t] + delij[t]) * inv_dists + Cdot_inv3 * (rsqik * delij[t] + rsqij * delik[t]);
    dcos_theta_dj[t] = delik[t] * inv_dists - Cdot_inv3 * rsqik * delij[t];
    dcos_theta_dk[t] = delij[t] * inv_dists - Cdot_inv3 * rsqij * delik[t];
  }

  sin_theta = sin(theta);
  if (sin_theta < 1.0e-5) sin_theta = 1.0e-5;
  p_val1 = paramsthbp(jtype,itype,ktype).p_val1;

  // ANGLE ENERGY

  p_val1 = paramsthbp(jtype,itype,ktype).p_val1;
  p_val2 = paramsthbp(jtype,itype,ktype).p_val2;
  p_val4 = paramsthbp(jtype,itype,ktype).p_val4;
  p_val7 = paramsthbp(jtype,itype,ktype).p_val7;
  theta_00 = paramsthbp(jtype,itype,ktype).theta_00;

  exp3ij = exp(-p_val3 * pow(BOA_ij, p_val4));
  f7_ij = 1.0 - exp3ij;
  Cf7ij = p_val3 * p_val4 * pow(BOA_ij, p_val4 - 1.0) * exp3ij;
  exp3jk = exp(-p_val3 * pow(BOA_ik, p_val4));
  f7_jk = 1.0 - exp3jk;
  Cf7jk = p_val3 * p_val4 * pow(BOA_ik, p_val4 - 1.0) * exp3jk;
  expval7 = exp(-p_val7 * d_Delta_boc[i]);
  trm8 = 1.0 + expval6 + expval7;
  f8_Dj = p_val5 - ((p_val5 - 1.0) * (2.0 + expval6) / trm8);
  Cf8j = ((1.0 - p_val5) / (trm8*trm8)) *
   (p_val6 * expval6 * trm8 - (2.0 + expval6) * (p_val6*expval6 - p_val7*expval7));
  theta_0 = 180.0 - theta_00 * (1.0 - exp(-p_val10 * (2.0 - SBO2)));
  theta_0 = theta_0*constPI/180.0;

  expval2theta  = exp(-p_val2 * (theta_0-theta)*(theta_0-theta));
  if (p_val1 >= 0)
    expval12theta = p_val1 * (1.0 - expval2theta);
  else // To avoid linear Me-H-Me angles (6/6/06)
    expval12theta = p_val1 * -expval2theta;

  CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
  CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
  CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
  CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj * expval2theta * (theta_0 - theta);
  Ctheta_0 = p_val10 * theta_00*constPI/180.0 * exp(-p_val10 * (2.0 - SBO2));
  CEval5 = -CEval4 * Ctheta_0 * CSBO2;
  CEval6 = CEval5 * dSBO1;
  CEval7 = CEval5 * dSBO2;
  CEval8 = -CEval4 / sin_theta;

  e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
  if (eflag) ev.ereax[3] += e_ang;

  // Penalty energy

  p_pen1 = paramsthbp(jtype,itype,ktype).p_pen1;

  exp_pen2ij = exp(-p_pen2 * (BOA_ij - 2.0)*(BOA_ij - 2.0));
  exp_pen2jk = exp(-p_pen2 * (BOA_ik - 2.0)*(BOA_ik - 2.0));
  exp_pen3 = exp(-p_pen3 * d_Delta[i]);
  exp_pen4 = exp(p_pen4 * d_Delta[i]);
  trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
  f9_Dj = (2.0 + exp_pen3) / trm_pen34;
  Cf9j = (-p_pen3 * exp_pen3 * trm_pen34 - (2.0 + exp_pen3) *
   (-p_pen3 * exp_pen3 + p_pen4 * exp_pen4))/(trm_pen34*trm_pen34);

  e_pen = p_pen1 * f9_Dj * exp_pen2ij * exp_pen2jk;
  if (eflag) ev.ereax[4] += e_pen;

  CEpen1 = e_pen * Cf9j / f9_Dj;
  temp   = -2.0 * p_pen2 * e_pen;
  CEpen2 = temp * (BOA_ij - 2.0);
  CEpen3 = temp * (BOA_ik - 2.0);

  // ConjAngle energy

  p_coa1 = paramsthbp(jtype,itype,ktype).p_coa1;
  exp_coa2 = exp(p_coa2 * Delta_val);
  e_coa = p_coa1 / (1. + exp_coa2) *
          exp(-p_coa3 * SQR(d_total_bo[j]-BOA_ij)) *
          exp(-p_coa3 * SQR(d_total_bo[k]-BOA_ik)) *
          exp(-p_coa4 * SQR(BOA_ij - 1.5)) *
          exp(-p_coa4 * SQR(BOA_ik - 1.5));

  CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
  CEcoa2 = -2 * p_coa4 * (BOA_ik - 1.5) * e_coa;
  CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
  CEcoa4 = -2 * p_coa3 * (d_total_bo[j]-BOA_ij) * e_coa;
  CEcoa5 = -2 * p_coa3 * (d_total_bo[k]-BOA_ik) * e_coa;

  if (eflag) ev.ereax[5] += e_coa;

  // Forces

  a_Cdbo(i,j_index) += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
  a_Cdbo(j,i_index) += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
  a_Cdbo(i,k_index) += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
  a_Cdbo(k,i_index) += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));

  CdDelta_i += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
  CdDelta_j += CEcoa4;
  a_CdDelta[k] += CEcoa5;

  for (int ll = j_start; ll < j_end; ll++) {
    int l = d_bo_list[ll];
    l &= NEIGHMASK;
    const int l_index = ll - j_start;

    temp_bo_jt = d_BO(i,l_index);
    temp = temp_bo_jt * temp_bo_jt * temp_bo_jt;
    pBOjt7 = temp * temp * temp_bo_jt;

    a_Cdbo(i,l_index) += (CEval6 * pBOjt7);
    a_Cdbopi(i,l_index) += CEval5;
    a_Cdbopi2(i,l_index) += CEval5;
  }

  for (int d = 0; d < 3; d++) fi_tmp[d] = CEval8 * dcos_theta_di[d];
  for (int d = 0; d < 3; d++) fj_tmp[d] = CEval8 * dcos_theta_dj[d];
  for (int d = 0; d < 3; d++) fk_tmp[d] = CEval8 * dcos_theta_dk[d];
  for (int d = 0; d < 3; d++) fitmp[d] -= fi_tmp[d];
  for (int d = 0; d < 3; d++) fjtmp[d] -= fj_tmp[d];
  for (int d = 0; d < 3; d++) a_f(k,d) -= fk_tmp[d];

  // energy/virial tally
  if (EVFLAG) {
    eng_tmp = e_ang + e_pen + e_coa;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng_tmp,0.0,0.0,0.0,0.0);
    for (int d = 0; d < 3; d++) delki[d] = -1.0 * delik[d];
    for (int d = 0; d < 3; d++) delji[d] = -1.0 * delij[d];
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
    if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj_tmp,fk_tmp,delji,delki);
  }

  a_CdDelta[j] += CdDelta_j;
  for (int d = 0; d < 3; d++) a_f(j,d) += fjtmp[d];
  a_CdDelta[i] += CdDelta_i;
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int &apack) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>(), apack, ev);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int &tpack, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  //auto a_Cdbo = dup_Cdbo.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // in reaxff_torsion_angles: j = i, k = j, i = k;

  F_FLOAT Delta_i, Delta_j, bo_ij, bo_ik, bo_jl, BOA_ij, BOA_ik, BOA_jl;
  F_FLOAT p_tor1, p_cot1, V1, V2, V3;
  F_FLOAT exp_tor2_ij, exp_tor2_ik, exp_tor2_jl, exp_tor1, exp_tor3_DiDj, exp_tor4_DiDj, exp_tor34_inv;
  F_FLOAT exp_cot2_ij, exp_cot2_ik, exp_cot2_jl, fn10, f11_DiDj, dfn11, fn12;
  F_FLOAT theta_ijk, theta_jil, sin_ijk, sin_jil, cos_ijk, cos_jil, tan_ijk_i, tan_jil_i;
  F_FLOAT cos_omega, cos2omega, cos3omega;
  F_FLOAT CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  F_FLOAT CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  F_FLOAT Cconj, CEconj1, CEconj2, CEconj3, CEconj4, CEconj5, CEconj6;
  F_FLOAT e_tor, e_con, eng_tmp;

  F_FLOAT delij[3], delik[3], deljl[3], dellk[3], delil[3], delkl[3];
  F_FLOAT fi_tmp[3], fj_tmp[3], fk_tmp[3];
  F_FLOAT dcos_ijk_di[3], dcos_ijk_dj[3], dcos_ijk_dk[3], dcos_jil_di[3], dcos_jil_dj[3], dcos_jil_dk[3];

  F_FLOAT p_tor2 = gp[23];
  F_FLOAT p_tor3 = gp[24];
  F_FLOAT p_tor4 = gp[25];
  F_FLOAT p_cot2 = gp[27];

  reax_int4 pack = d_torsion_pack(tpack,0);
  const int i = pack.i0;
  const int j = pack.i1;
  const int k = pack.i2;
  const int l = pack.i3;

  pack = d_torsion_pack(tpack, 1);
  //const int i = pack.i0;
  const int j_index = pack.i1;
  const int k_index = pack.i2;
  const int l_index = pack.i3;

  const int itype = type(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  Delta_i = d_Delta_boc[i];

  const int jtype = type(j);

  bo_ij = d_BO(i,j_index);

  delij[0] = x(j,0) - xtmp;
  delij[1] = x(j,1) - ytmp;
  delij[2] = x(j,2) - ztmp;
  const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
  const F_FLOAT rij = sqrt(rsqij);

  BOA_ij = bo_ij - thb_cut;
  Delta_j = d_Delta_boc[j];
  exp_tor2_ij = exp(-p_tor2 * BOA_ij);
  exp_cot2_ij = exp(-p_cot2 * SQR(BOA_ij - 1.5));
  exp_tor3_DiDj = exp(-p_tor3 * (Delta_i + Delta_j));
  exp_tor4_DiDj = exp(p_tor4  * (Delta_i + Delta_j));
  exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DiDj + exp_tor4_DiDj);
  f11_DiDj = (2.0 + exp_tor3_DiDj) * exp_tor34_inv;

  const int ktype = type(k);

  bo_ik = d_BO(i,k_index);

  BOA_ik = bo_ik - thb_cut;
  for (int d = 0; d < 3; d ++) delik[d] = x(k,d) - x(i,d);
  const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
  const F_FLOAT rik = sqrt(rsqik);

  cos_ijk = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
  if (cos_ijk > 1.0) cos_ijk  = 1.0;
  if (cos_ijk < -1.0) cos_ijk  = -1.0;
  theta_ijk = acos(cos_ijk);

  // dcos_ijk
  const F_FLOAT inv_dists = 1.0 / (rij * rik);
  const F_FLOAT cos_ijk_tmp = cos_ijk / ((rij*rik)*(rij*rik));

  for (int d = 0; d < 3; d++) {
    dcos_ijk_di[d] = -(delik[d] + delij[d]) * inv_dists + cos_ijk_tmp * (rsqik * delij[d] + rsqij * delik[d]);
    dcos_ijk_dj[d] = delik[d] * inv_dists - cos_ijk_tmp * rsqik * delij[d];
    dcos_ijk_dk[d] = delij[d] * inv_dists - cos_ijk_tmp * rsqij * delik[d];
  }

  sin_ijk = sin(theta_ijk);
  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE)
    tan_ijk_i = cos_ijk / MIN_SINE;
  else if (sin_ijk <= 0 && sin_ijk >= -MIN_SINE)
    tan_ijk_i = -cos_ijk / MIN_SINE;
  else tan_ijk_i = cos_ijk / sin_ijk;

  exp_tor2_ik = exp(-p_tor2 * BOA_ik);
  exp_cot2_ik = exp(-p_cot2 * SQR(BOA_ik -1.5));

  const int ltype = type(l);

  bo_jl = d_BO(j,l_index);

  for (int d = 0; d < 3; d ++) deljl[d] = x(l,d) - x(j,d);
  const F_FLOAT rsqjl = deljl[0]*deljl[0] + deljl[1]*deljl[1] + deljl[2]*deljl[2];
  const F_FLOAT rjl = sqrt(rsqjl);
  BOA_jl = bo_jl - thb_cut;

  cos_jil = -(delij[0]*deljl[0]+delij[1]*deljl[1]+delij[2]*deljl[2])/(rij*rjl);
  if (cos_jil > 1.0) cos_jil  = 1.0;
  if (cos_jil < -1.0) cos_jil  = -1.0;
  theta_jil = acos(cos_jil);

  // dcos_jil
  const F_FLOAT inv_distjl = 1.0 / (rij * rjl);
  const F_FLOAT cos_jil_tmp = cos_jil / ((rij*rjl)*(rij*rjl));

  for (int d = 0; d < 3; d++) {
    dcos_jil_di[d] = deljl[d] * inv_distjl - cos_jil_tmp * rsqjl * -delij[d];
    dcos_jil_dj[d] = (-deljl[d] + delij[d]) * inv_distjl - cos_jil_tmp * (rsqjl * delij[d] + rsqij * -deljl[d]);
    dcos_jil_dk[d] = -delij[d] * inv_distjl - cos_jil_tmp * rsqij * deljl[d];
  }

  sin_jil = sin(theta_jil);
  if (sin_jil >= 0 && sin_jil <= MIN_SINE)
    tan_jil_i = cos_jil / MIN_SINE;
  else if (sin_jil <= 0 && sin_jil >= -MIN_SINE)
    tan_jil_i = -cos_jil / MIN_SINE;
  else tan_jil_i = cos_jil / sin_jil;

  for (int d = 0; d < 3; d ++) dellk[d] = x(k,d) - x(l,d);
  const F_FLOAT rsqlk = dellk[0]*dellk[0] + dellk[1]*dellk[1] + dellk[2]*dellk[2];
  const F_FLOAT rlk = sqrt(rsqlk);

  // non-Kokkos ReaxFF has a separate function for computing omega, which
  //  limits the scope of the MIN_SINE statements below

  F_FLOAT sin_ijk_rnd = sin_ijk;
  F_FLOAT sin_jil_rnd = sin_jil;

  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE) sin_ijk_rnd = MIN_SINE;
  else if (sin_ijk <= 0 && sin_ijk >= -MIN_SINE) sin_ijk_rnd = -MIN_SINE;
  if (sin_jil >= 0 && sin_jil <= MIN_SINE) sin_jil_rnd = MIN_SINE;
  else if (sin_jil <= 0 && sin_jil >= -MIN_SINE) sin_jil_rnd = -MIN_SINE;

  F_FLOAT unnorm_cos_omega, unnorm_sin_omega, omega;
  F_FLOAT htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  F_FLOAT arg, poem, tel;
  F_FLOAT cross_ij_jl[3];

  // omega

  F_FLOAT dot_ij_jk = -(delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2]);
  F_FLOAT dot_ij_lj = delij[0]*deljl[0]+delij[1]*deljl[1]+delij[2]*deljl[2];
  F_FLOAT dot_ik_jl = delik[0]*deljl[0]+delik[1]*deljl[1]+delik[2]*deljl[2];
  unnorm_cos_omega = dot_ij_jk * dot_ij_lj + rsqij * dot_ik_jl;

  cross_ij_jl[0] = delij[1]*deljl[2] - delij[2]*deljl[1];
  cross_ij_jl[1] = delij[2]*deljl[0] - delij[0]*deljl[2];
  cross_ij_jl[2] = delij[0]*deljl[1] - delij[1]*deljl[0];

  unnorm_sin_omega = -rij*(delik[0]*cross_ij_jl[0]+delik[1]*cross_ij_jl[1]+delik[2]*cross_ij_jl[2]);
  omega = atan2(unnorm_sin_omega, unnorm_cos_omega);

  htra = rik + cos_ijk * (rjl * cos_jil - rij);
  htrb = rij - rik * cos_ijk - rjl * cos_jil;
  htrc = rjl + cos_jil * (rik * cos_ijk - rij);
  hthd = rik * sin_ijk_rnd * (rij - rjl * cos_jil);
  hthe = rjl * sin_jil_rnd * (rij - rik * cos_ijk);
  hnra = rjl * sin_ijk_rnd * sin_jil_rnd;
  hnrc = rik * sin_ijk_rnd * sin_jil_rnd;
  hnhd = rik * rjl * cos_ijk * sin_jil_rnd;
  hnhe = rik * rjl * sin_ijk_rnd * cos_jil;

  tel = SQR(rik) + SQR(rij) + SQR(rjl) - SQR(rlk) -
        2.0 * (rik * rij * cos_ijk - rik * rjl * cos_ijk * cos_jil + rij * rjl * cos_jil);

  poem = 2.0 * rik * rjl * sin_ijk_rnd * sin_jil_rnd;
  F_FLOAT inv_poem = 1.0 / poem;

  arg = tel * inv_poem;
  if (arg >  1.0) arg =  1.0;
  if (arg < -1.0) arg = -1.0;

  cos_omega = cos(omega);
  cos2omega = cos(2. * omega);
  cos3omega = cos(3. * omega);

  // torsion energy

  p_tor1 = paramsfbp(ktype,itype,jtype,ltype).p_tor1;
  p_cot1 = paramsfbp(ktype,itype,jtype,ltype).p_cot1;
  V1 = paramsfbp(ktype,itype,jtype,ltype).V1;
  V2 = paramsfbp(ktype,itype,jtype,ltype).V2;
  V3 = paramsfbp(ktype,itype,jtype,ltype).V3;

  exp_tor1 = exp(p_tor1 * SQR(2.0 - d_BO_pi(i,j_index) - f11_DiDj));
  exp_tor2_jl = exp(-p_tor2 * BOA_jl);
  exp_cot2_jl = exp(-p_cot2 * SQR(BOA_jl - 1.5));
  fn10 = (1.0 - exp_tor2_ik) * (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jl);

  CV = 0.5 * (V1 * (1.0 + cos_omega) + V2 * exp_tor1 * (1.0 - cos2omega) + V3 * (1.0 + cos3omega));

  e_tor = fn10 * sin_ijk * sin_jil * CV;
  if (eflag) ev.ereax[6] += e_tor;

  dfn11 = (-p_tor3 * exp_tor3_DiDj + (p_tor3 * exp_tor3_DiDj - p_tor4 * exp_tor4_DiDj) *
          (2.0 + exp_tor3_DiDj) * exp_tor34_inv) * exp_tor34_inv;

  CEtors1 = sin_ijk * sin_jil * CV;

  CEtors2 = -fn10 * 2.0 * p_tor1 * V2 * exp_tor1 * (2.0 - d_BO_pi(i,j_index) - f11_DiDj) *
            (1.0 - SQR(cos_omega)) * sin_ijk * sin_jil;
  CEtors3 = CEtors2 * dfn11;

  CEtors4 = CEtors1 * p_tor2 * exp_tor2_ik * (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jl);
  CEtors5 = CEtors1 * p_tor2 * (1.0 - exp_tor2_ik) * exp_tor2_ij * (1.0 - exp_tor2_jl);
  CEtors6 = CEtors1 * p_tor2 * (1.0 - exp_tor2_ik) * (1.0 - exp_tor2_ij) * exp_tor2_jl;

  cmn = -fn10 * CV;
  CEtors7 = cmn * sin_jil * tan_ijk_i;
  CEtors8 = cmn * sin_ijk * tan_jil_i;

  CEtors9 = fn10 * sin_ijk * sin_jil *
    (0.5 * V1 - 2.0 * V2 * exp_tor1 * cos_omega + 1.5 * V3 * (cos2omega + 2.0 * SQR(cos_omega)));

  // 4-body conjugation energy

  fn12 = exp_cot2_ik * exp_cot2_ij * exp_cot2_jl;
  e_con = p_cot1 * fn12 * (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jil);
  if (eflag) ev.ereax[7] += e_con;

  Cconj = -2.0 * fn12 * p_cot1 * p_cot2 * (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jil);

  CEconj1 = Cconj * (BOA_ik - 1.5e0);
  CEconj2 = Cconj * (BOA_ij - 1.5e0);
  CEconj3 = Cconj * (BOA_jl - 1.5e0);

  CEconj4 = -p_cot1 * fn12 * (SQR(cos_omega) - 1.0) * sin_jil * tan_ijk_i;
  CEconj5 = -p_cot1 * fn12 * (SQR(cos_omega) - 1.0) * sin_ijk * tan_jil_i;
  CEconj6 = 2.0 * p_cot1 * fn12 * cos_omega * sin_ijk * sin_jil;

  // forces

  // contribution to bond order

  a_Cdbopi(i,j_index) += CEtors2;

  a_CdDelta[j] += CEtors3;
  a_CdDelta[i] += CEtors3;

  a_Cdbo(i,k_index) += CEtors4 + CEconj1;
  a_Cdbo(i,j_index) += CEtors5 + CEconj2;
  a_Cdbo(j,l_index) += CEtors6 + CEconj3;

  const F_FLOAT coeff74 = CEtors7 + CEconj4;
  const F_FLOAT coeff85 = CEtors8 + CEconj5;
  const F_FLOAT coeff96 = CEtors9 + CEconj6;

  const F_FLOAT inv_rij = 1.0 / rij;
  const F_FLOAT inv_rik = 1.0 / rik;
  const F_FLOAT inv_rjl = 1.0 / rjl;
  const F_FLOAT inv_sin_ijk_rnd = 1.0 / sin_ijk_rnd;
  const F_FLOAT inv_sin_jil_rnd = 1.0 / sin_jil_rnd;

#ifdef LMP_KOKKOS_GPU
  #pragma unroll
#endif
  for (int d = 0; d < 3; d++) {
    // dcos_omega_di
    F_FLOAT dcos_omega_dk = ((htra-arg*hnra) * inv_rik) * delik[d] - dellk[d];
    dcos_omega_dk += (hthd-arg*hnhd) * inv_sin_ijk_rnd * -dcos_ijk_dk[d];
    dcos_omega_dk *= 2.0 * inv_poem;

    // dcos_omega_dj
    F_FLOAT dcos_omega_di = -((htra-arg*hnra) * inv_rik) * delik[d] - htrb * inv_rij * delij[d];
    dcos_omega_di += -(hthd-arg*hnhd) * inv_sin_ijk_rnd * dcos_ijk_di[d];
    dcos_omega_di += -(hthe-arg*hnhe) * inv_sin_jil_rnd * dcos_jil_di[d];
    dcos_omega_di *= 2.0 * inv_poem;

    // dcos_omega_dk
    F_FLOAT dcos_omega_dj = -((htrc-arg*hnrc) * inv_rjl) * deljl[d] + htrb * inv_rij * delij[d];
    dcos_omega_dj += -(hthd-arg*hnhd) * inv_sin_ijk_rnd * dcos_ijk_dj[d];
    dcos_omega_dj += -(hthe-arg*hnhe) * inv_sin_jil_rnd * dcos_jil_dj[d];
    dcos_omega_dj *= 2.0 * inv_poem;

    // dcos_omega_dl
    F_FLOAT dcos_omega_dl = ((htrc-arg*hnrc) * inv_rjl) * deljl[d] + dellk[d];
    dcos_omega_dl += (hthe-arg*hnhe) * inv_sin_jil_rnd * -dcos_jil_dk[d];
    dcos_omega_dl *= 2.0 * inv_poem;

    // dcos_theta_ijk
    fi_tmp[d] = (coeff74) * dcos_ijk_di[d];
    fj_tmp[d] = (coeff74) * dcos_ijk_dj[d];
    fk_tmp[d] = (coeff74) * dcos_ijk_dk[d];

    // dcos_theta_jil
    fi_tmp[d] += (coeff85) * dcos_jil_di[d];
    fj_tmp[d] += (coeff85) * dcos_jil_dj[d];
    F_FLOAT fl_tmp =  (coeff85) * dcos_jil_dk[d];

    // dcos_omega
    fi_tmp[d] += (coeff96) * dcos_omega_di;
    fj_tmp[d] += (coeff96) * dcos_omega_dj;
    fk_tmp[d] += (coeff96) * dcos_omega_dk;
    fl_tmp += (coeff96) * dcos_omega_dl;

    // total forces
    a_f(i,d) -= fi_tmp[d];
    a_f(j,d) -= fj_tmp[d];
    a_f(k,d) -= fk_tmp[d];
    a_f(l,d) -= fl_tmp;
  }

  // per-atom energy/virial tally

  if (EVFLAG) {
    eng_tmp = e_tor + e_con;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng_tmp,0.0,0.0,0.0,0.0);
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
    if (vflag_either) {
        for (int d = 0; d < 3; d ++) delil[d] = x(l,d) - x(i,d);
        for (int d = 0; d < 3; d ++) delkl[d] = x(l,d) - x(k,d);
        this->template v_tally4<NEIGHFLAG>(ev,k,i,j,l,fk_tmp,fi_tmp,fj_tmp,delkl,delil,deljl);
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int &tpack) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>(), tpack, ev);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int hblist[MAX_BONDS];
  F_FLOAT theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
  F_FLOAT e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
  F_FLOAT dcos_theta_di[3], dcos_theta_dj[3], dcos_theta_dk[3];

  // tally variables
  F_FLOAT fi_tmp[3], fj_tmp[3], fk_tmp[3], delij[3], delji[3], delik[3], delki[3];
  for (int d = 0; d < 3; d++) fi_tmp[d] = fj_tmp[d] = fk_tmp[d] = 0.0;

  const int i = d_ilist[ii];
  const int itype = type(i);
  if (paramssing(itype).p_hbond != 1) return;

  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];
  const int k_start = d_hb_first[i];
  const int k_end = k_start + d_hb_num[i];

  int top = 0;
  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const int jtype = type(j);
    const int j_index = jj - j_start;
    const F_FLOAT bo_ij = d_BO(i,j_index);

    if (paramssing(jtype).p_hbond == 2 && bo_ij >= HB_THRESHOLD) {
      hblist[top] = jj;
      top ++;
    }
  }

  F_FLOAT fitmp[3];
  for (int d = 0; d < 3; d++) fitmp[d] = 0.0;

  for (int kk = k_start; kk < k_end; kk++) {
    int k = d_hb_list[kk];
    k &= NEIGHMASK;
    const tagint ktag = tag(k);
    const int ktype = type(k);

    delik[0] = x(k,0) - xtmp;
    delik[1] = x(k,1) - ytmp;
    delik[2] = x(k,2) - ztmp;
    const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
    const F_FLOAT rik = sqrt(rsqik);

    for (int itr = 0; itr < top; itr++) {
      const int jj = hblist[itr];
      int j = d_bo_list[jj];
      j &= NEIGHMASK;
      const tagint jtag = tag(j);
      if (jtag == ktag) continue;

      const int jtype = type(j);
      const int j_index = jj - j_start;
      const F_FLOAT bo_ij = d_BO(i,j_index);

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
      const F_FLOAT rij = sqrt(rsqij);

      // theta and derivatives
      cos_theta = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
      if (cos_theta > 1.0) cos_theta  = 1.0;
      if (cos_theta < -1.0) cos_theta  = -1.0;
      theta = acos(cos_theta);

      const F_FLOAT inv_dists = 1.0 / (rij * rik);
      const F_FLOAT Cdot_inv3 = cos_theta * inv_dists * inv_dists;

      for (int d = 0; d < 3; d++) {
        dcos_theta_di[d] = -(delik[d] + delij[d]) * inv_dists + Cdot_inv3 * (rsqik * delij[d] + rsqij * delik[d]);
        dcos_theta_dj[d] = delik[d] * inv_dists - Cdot_inv3 * rsqik * delij[d];
        dcos_theta_dk[d] = delij[d] * inv_dists - Cdot_inv3 * rsqij * delik[d];
      }

      // hydrogen bond energy
      const F_FLOAT p_hb1 = paramshbp(jtype,itype,ktype).p_hb1;
      const F_FLOAT p_hb2 = paramshbp(jtype,itype,ktype).p_hb2;
      const F_FLOAT p_hb3 = paramshbp(jtype,itype,ktype).p_hb3;
      const F_FLOAT r0_hb = paramshbp(jtype,itype,ktype).r0_hb;

      sin_theta2 = sin(theta/2.0);
      sin_xhz4 = SQR(sin_theta2);
      sin_xhz4 *= sin_xhz4;
      cos_xhz1 = (1.0 - cos_theta);
      exp_hb2 = exp(-p_hb2 * bo_ij);
      exp_hb3 = exp(-p_hb3 * (r0_hb/rik + rik/r0_hb - 2.0));

      e_hb = p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;
      if (EVFLAG && eflag_global) ev.ereax[8] += e_hb;

      // hydrogen bond forces
      CEhb1 = p_hb1 * p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
      CEhb2 = -p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
      CEhb3 = -p_hb3 * (-r0_hb/SQR(rik) + 1.0/r0_hb) * e_hb;

      d_Cdbo(i,j_index) += CEhb1; // dbo term

      // dcos terms
      for (int d = 0; d < 3; d++) fi_tmp[d] = CEhb2 * dcos_theta_di[d];
      for (int d = 0; d < 3; d++) fj_tmp[d] = CEhb2 * dcos_theta_dj[d];
      for (int d = 0; d < 3; d++) fk_tmp[d] = CEhb2 * dcos_theta_dk[d];

      // dr terms
      for (int d = 0; d < 3; d++) fi_tmp[d] -= CEhb3/rik * delik[d];
      for (int d = 0; d < 3; d++) fk_tmp[d] += CEhb3/rik * delik[d];

      for (int d = 0; d < 3; d++) fitmp[d] -= fi_tmp[d];
      for (int d = 0; d < 3; d++) a_f(j,d) -= fj_tmp[d];
      for (int d = 0; d < 3; d++) a_f(k,d) -= fk_tmp[d];

      for (int d = 0; d < 3; d++) delki[d] = -1.0 * delik[d];
      for (int d = 0; d < 3; d++) delji[d] = -1.0 * delij[d];

      if (EVFLAG) {
        if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,e_hb);
        if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj_tmp,fk_tmp,delji,delki);
      }
    }
  }
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>(), ii, ev);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxUpdateBond<NEIGHFLAG>, const int &ii) const {

  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  Kokkos::View<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi2 = d_Cdbopi2;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const tagint itag = tag(i);
  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    int flag = 0;

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) flag = 1;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) flag = 1;
    } else {
      if (x(j,2)  < ztmp) flag = 1;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) flag = 1;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) flag = 1;
    }

    if (!flag) continue;

    const int j_index = jj - j_start;
    const F_FLOAT Cdbo_i = d_Cdbo(i,j_index);
    const F_FLOAT Cdbopi_i = d_Cdbopi(i,j_index);
    const F_FLOAT Cdbopi2_i = d_Cdbopi2(i,j_index);

    const int k_start = d_bo_first[j];
    const int k_end = k_start + d_bo_num[j];

    for (int kk = k_start; kk < k_end; kk++) {
      int k = d_bo_list[kk];
      k &= NEIGHMASK;
      if (k != i) continue;
      const int k_index = kk - k_start;

      a_Cdbo(j,k_index) += Cdbo_i;
      a_Cdbopi(j,k_index) += Cdbopi_i;
      a_Cdbopi2(j,k_index) += Cdbopi2_i;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT p_be1, p_be2, De_s, De_p, De_pp, pow_BOs_be2, exp_be12, CEbo, ebond;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const tagint itag = tag(i);
  const F_FLOAT imass = paramssing(itype).mass;
  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  F_FLOAT CdDelta_i = 0.0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const int jtype = type(j);
    const int j_index = jj - j_start;
    const F_FLOAT jmass = paramssing(jtype).mass;

    // bond energy (nlocal only)
    p_be1 = paramstwbp(itype,jtype).p_be1;
    p_be2 = paramstwbp(itype,jtype).p_be2;
    De_s = paramstwbp(itype,jtype).De_s;
    De_p = paramstwbp(itype,jtype).De_p;
    De_pp = paramstwbp(itype,jtype).De_pp;

    const F_FLOAT BO_i = d_BO(i,j_index);
    const F_FLOAT BO_s_i = d_BO_s(i,j_index);
    const F_FLOAT BO_pi_i = d_BO_pi(i,j_index);
    const F_FLOAT BO_pi2_i = d_BO_pi2(i,j_index);

    if (BO_s_i == 0.0) pow_BOs_be2 = 0.0;
    else pow_BOs_be2 = pow(BO_s_i,p_be2);
    exp_be12 = exp(p_be1*(1.0-pow_BOs_be2));
    CEbo = -De_s*exp_be12*(1.0-p_be1*p_be2*pow_BOs_be2);
    ebond = -De_s*BO_s_i*exp_be12
                              -De_p*BO_pi_i
                          -De_pp*BO_pi2_i;

    if (EFLAG && eflag_global) ev.evdwl += ebond;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,ebond,0.0,0.0,0.0,0.0);
    //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,ebond);

    // calculate derivatives of Bond Orders
    d_Cdbo(i,j_index) += CEbo;
    d_Cdbopi(i,j_index) -= (CEbo + De_p);
    d_Cdbopi2(i,j_index) -= (CEbo + De_pp);

    // Stabilisation terminal triple bond
    F_FLOAT estriph = 0.0;

    if (BO_i >= 1.00) {
      if (gp[37] == 2 || (imass == 12.0000 && jmass == 15.9990) ||
                         (jmass == 12.0000 && imass == 15.9990)) {
        const F_FLOAT exphu = exp(-gp[7] * SQR(BO_i - 2.50));
        const F_FLOAT exphua1 = exp(-gp[3] * (d_total_bo[i]-BO_i));
        const F_FLOAT exphub1 = exp(-gp[3] * (d_total_bo[j]-BO_i));
        const F_FLOAT exphuov = exp(gp[4] * (d_Delta[i] + d_Delta[j]));
        const F_FLOAT hulpov = 1.0 / (1.0 + 25.0 * exphuov);
        estriph = gp[10] * exphu * hulpov * (exphua1 + exphub1);

        if (EFLAG && eflag_global) ev.evdwl += estriph;
        //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,estriph,0.0,0.0,0.0,0.0);
        //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,estriph);

        const F_FLOAT decobdbo = gp[10] * exphu * hulpov * (exphua1 + exphub1) *
            (gp[3] - 2.0 * gp[7] * (BO_i-2.50));
        const F_FLOAT decobdboua = -gp[10] * exphu * hulpov *
            (gp[3]*exphua1 + 25.0*gp[4]*exphuov*hulpov*(exphua1+exphub1));
        const F_FLOAT decobdboub = -gp[10] * exphu * hulpov *
            (gp[3]*exphub1 + 25.0*gp[4]*exphuov*hulpov*(exphua1+exphub1));

        d_Cdbo(i,j_index) += decobdbo;
        CdDelta_i += decobdboua;
        a_CdDelta[j] += decobdboub;
      }
    }
    const F_FLOAT eng_tmp = ebond + estriph;
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
  }
  a_CdDelta[i] += CdDelta_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EFLAG>(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int VFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT delij[3], delik[3], deljk[3], tmpvec[3];
  F_FLOAT dBOp_i[3], dBOp_k[3], dln_BOp_pi[3], dln_BOp_pi2[3];

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const tagint itag = tag(i);
  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];

  F_FLOAT CdDelta_i = d_CdDelta[i];
  F_FLOAT fitmp[3];
  for (int j = 0; j < 3; j++) fitmp[j] = 0.0;

  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const int j_index = jj - j_start;
    F_FLOAT CdDelta_j = d_CdDelta[j];

    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;

    const int k_start = d_bo_first[j];
    const int k_end = k_start + d_bo_num[j];

    F_FLOAT coef_C1dbo, coef_C2dbo, coef_C3dbo, coef_C1dbopi, coef_C2dbopi, coef_C3dbopi, coef_C4dbopi;
    F_FLOAT coef_C1dbopi2, coef_C2dbopi2, coef_C3dbopi2, coef_C4dbopi2, coef_C1dDelta, coef_C2dDelta, coef_C3dDelta;

    coef_C1dbo = coef_C2dbo = coef_C3dbo = 0.0;
    coef_C1dbopi = coef_C2dbopi = coef_C3dbopi = coef_C4dbopi = 0.0;
    coef_C1dbopi2 = coef_C2dbopi2 = coef_C3dbopi2 = coef_C4dbopi2 = 0.0;
    coef_C1dDelta = coef_C2dDelta = coef_C3dDelta = 0.0;

    // total forces on i, j, k (nlocal + nghost, from Add_dBond_to_Forces))
    const F_FLOAT Cdbo_ij = d_Cdbo(i,j_index);
    coef_C1dbo = d_C1dbo(i,j_index) * (Cdbo_ij);
    coef_C2dbo = d_C2dbo(i,j_index) * (Cdbo_ij);
    coef_C3dbo = d_C3dbo(i,j_index) * (Cdbo_ij);

    const F_FLOAT Cdbopi_ij = d_Cdbopi(i,j_index);
    coef_C1dbopi = d_C1dbopi(i,j_index) * (Cdbopi_ij);
    coef_C2dbopi = d_C2dbopi(i,j_index) * (Cdbopi_ij);
    coef_C3dbopi = d_C3dbopi(i,j_index) * (Cdbopi_ij);
    coef_C4dbopi = d_C4dbopi(i,j_index) * (Cdbopi_ij);

    const F_FLOAT Cdbopi2_ij = d_Cdbopi2(i,j_index);
    coef_C1dbopi2 = d_C1dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C2dbopi2 = d_C2dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C3dbopi2 = d_C3dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C4dbopi2 = d_C4dbopi2(i,j_index) * (Cdbopi2_ij);

    const F_FLOAT coeff_CdDelta_ij = CdDelta_i + CdDelta_j;
    coef_C1dDelta = d_C1dbo(i,j_index) * (coeff_CdDelta_ij);
    coef_C2dDelta = d_C2dbo(i,j_index) * (coeff_CdDelta_ij);
    coef_C3dDelta = d_C3dbo(i,j_index) * (coeff_CdDelta_ij);

    F_FLOAT temp[3];

    F_FLOAT d_dln_BOp_pi_local = d_dln_BOp_pi(i,j_index);
    for (int d = 0; d < 3; d++) dln_BOp_pi[d] = d_dln_BOp_pi_local * delij[d];

    F_FLOAT d_dln_BOp_pi2_local = d_dln_BOp_pi2(i,j_index);
    for (int d = 0; d < 3; d++) dln_BOp_pi2[d] = d_dln_BOp_pi2_local * delij[d];

    F_FLOAT d_dBOp_local = d_dBOp(i,j_index);
    for (int d = 0; d < 3; d++) dBOp_i[d] = d_dBOp_local * delij[d];

    // forces on i
    for (int d = 0; d < 3; d++) temp[d] =  coef_C1dbo * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbo * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dDelta * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dDelta * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dbopi * dln_BOp_pi[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbopi * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbopi * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dbopi2 * dln_BOp_pi2[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbopi2 * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbopi2 * d_dDeltap_self(i,d);

    if (VFLAG && vflag_either) this->template v_tally<NEIGHFLAG>(ev,i,temp,delij);

    fitmp[0] -= temp[0];
    fitmp[1] -= temp[1];
    fitmp[2] -= temp[2];

    // forces on j
    for (int d = 0; d < 3; d++) temp[d] = -coef_C1dbo * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbo * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dDelta * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dDelta * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dbopi * dln_BOp_pi[d];
    for (int d = 0; d < 3; d++) temp[d] -= coef_C2dbopi * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C4dbopi * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dbopi2 * dln_BOp_pi2[d];
    for (int d = 0; d < 3; d++) temp[d] -= coef_C2dbopi2 * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C4dbopi2 * d_dDeltap_self(j,d);

    a_f(j,0) -= temp[0];
    a_f(j,1) -= temp[1];
    a_f(j,2) -= temp[2];

    if (VFLAG && vflag_either) {
      for (int d = 0; d < 3; d++) tmpvec[d] = -delij[d];
      this->template v_tally<NEIGHFLAG>(ev,j,temp,tmpvec);
    }

    // forces on k: i neighbor
    for (int kk = j_start; kk < j_end; kk++) {
      int k = d_bo_list[kk];
      k &= NEIGHMASK;
      const int k_index = kk - j_start;

      delik[0] = x(k,0) - xtmp;
      delik[1] = x(k,1) - ytmp;
      delik[2] = x(k,2) - ztmp;

      d_dBOp_local = d_dBOp(i,k_index);
      for (int d = 0; d < 3; d++) dBOp_k[d] = d_dBOp_local * delik[d];

      const F_FLOAT coef_all = -coef_C2dbo - coef_C2dDelta - coef_C3dbopi - coef_C3dbopi2;
      for (int d = 0; d < 3; d++) temp[d] = coef_all * dBOp_k[d];

      a_f(k,0) -= temp[0];
      a_f(k,1) -= temp[1];
      a_f(k,2) -= temp[2];

      if (VFLAG && vflag_either) {
        delik[0] = x(k,0) - xtmp;
        delik[1] = x(k,1) - ytmp;
        delik[2] = x(k,2) - ztmp;
        for (int d = 0; d < 3; d++) tmpvec[d] = x(j,d) - x(k,d) - delik[d];
        this->template v_tally<NEIGHFLAG>(ev,k,temp,tmpvec);
      }

    }

    // forces on k: j neighbor
    for (int kk = k_start; kk < k_end; kk++) {
      int k = d_bo_list[kk];
      k &= NEIGHMASK;
      const int k_index = kk - k_start;

      for (int d = 0; d < 3; d++) deljk[d] = x(k,d) - x(j,d);

      d_dBOp_local = d_dBOp(j,k_index);
      for (int d = 0; d < 3; d++) dBOp_k[d] = d_dBOp_local * deljk[d];
      const F_FLOAT coef_all = -coef_C3dbo - coef_C3dDelta - coef_C4dbopi - coef_C4dbopi2;
      for (int d = 0; d < 3; d++) temp[d] = coef_all * dBOp_k[d];

      a_f(k,0) -= temp[0];
      a_f(k,1) -= temp[1];
      a_f(k,2) -= temp[2];

      if (VFLAG && vflag_either) {
        for (int d = 0; d < 3; d++) deljk[d] = x(k,d) - x(j,d);
        for (int d = 0; d < 3; d++) tmpvec[d] = x(i,d) - x(k,d) - deljk[d];
        this->template v_tally<NEIGHFLAG>(ev,k,temp,tmpvec);
      }
    }
  }
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int VFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,VFLAG>(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::ev_tally(EV_FLOAT_REAX &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    a_eatom[i] += epairhalf;
    a_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      ev.v[0] += v0;
      ev.v[1] += v1;
      ev.v[2] += v2;
      ev.v[3] += v3;
      ev.v[4] += v4;
      ev.v[5] += v5;
    }

    if (vflag_atom) {
      a_vatom(i,0) += 0.5*v0;
      a_vatom(i,1) += 0.5*v1;
      a_vatom(i,2) += 0.5*v2;
      a_vatom(i,3) += 0.5*v3;
      a_vatom(i,4) += 0.5*v4;
      a_vatom(i,5) += 0.5*v5;
      a_vatom(j,0) += 0.5*v0;
      a_vatom(j,1) += 0.5*v1;
      a_vatom(j,2) += 0.5*v2;
      a_vatom(j,3) += 0.5*v3;
      a_vatom(j,4) += 0.5*v4;
      a_vatom(j,5) += 0.5*v5;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::e_tally(EV_FLOAT_REAX & /*ev*/, const int &i, const int &j,
      const F_FLOAT &epair) const
{
  // The eatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const E_FLOAT epairhalf = 0.5 * epair;
  a_eatom[i] += epairhalf;
  a_eatom[j] += epairhalf;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::e_tally_single(EV_FLOAT_REAX & /*ev*/, const int &i,
      const F_FLOAT &epair) const
{
  // The eatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  a_eatom[i] += epair;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally(EV_FLOAT_REAX &ev, const int &i,
  F_FLOAT *fi, F_FLOAT *drij) const
{
  F_FLOAT v[6];

  v[0] = 0.5*drij[0]*fi[0];
  v[1] = 0.5*drij[1]*fi[1];
  v[2] = 0.5*drij[2]*fi[2];
  v[3] = 0.5*drij[0]*fi[1];
  v[4] = 0.5*drij[0]*fi[2];
  v[5] = 0.5*drij[1]*fi[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    a_vatom(i,0) += v[0]; a_vatom(i,1) += v[1]; a_vatom(i,2) += v[2];
    a_vatom(i,3) += v[3]; a_vatom(i,4) += v[4]; a_vatom(i,5) += v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally3(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
  F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT v[6];

  v[0] = drij[0]*fj[0] + drik[0]*fk[0];
  v[1] = drij[1]*fj[1] + drik[1]*fk[1];
  v[2] = drij[2]*fj[2] + drik[2]*fk[2];
  v[3] = drij[0]*fj[1] + drik[0]*fk[1];
  v[4] = drij[0]*fj[2] + drik[0]*fk[2];
  v[5] = drij[1]*fj[2] + drik[1]*fk[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    a_vatom(i,0) += THIRD * v[0]; a_vatom(i,1) += THIRD * v[1]; a_vatom(i,2) += THIRD * v[2];
    a_vatom(i,3) += THIRD * v[3]; a_vatom(i,4) += THIRD * v[4]; a_vatom(i,5) += THIRD * v[5];
    a_vatom(j,0) += THIRD * v[0]; a_vatom(j,1) += THIRD * v[1]; a_vatom(j,2) += THIRD * v[2];
    a_vatom(j,3) += THIRD * v[3]; a_vatom(j,4) += THIRD * v[4]; a_vatom(j,5) += THIRD * v[5];
    a_vatom(k,0) += THIRD * v[0]; a_vatom(k,1) += THIRD * v[1]; a_vatom(k,2) += THIRD * v[2];
    a_vatom(k,3) += THIRD * v[3]; a_vatom(k,4) += THIRD * v[4]; a_vatom(k,5) += THIRD * v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally4(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
  const int &l, F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *dril, F_FLOAT *drjl, F_FLOAT *drkl) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  F_FLOAT v[6];

  v[0] = dril[0]*fi[0] + drjl[0]*fj[0] + drkl[0]*fk[0];
  v[1] = dril[1]*fi[1] + drjl[1]*fj[1] + drkl[1]*fk[1];
  v[2] = dril[2]*fi[2] + drjl[2]*fj[2] + drkl[2]*fk[2];
  v[3] = dril[0]*fi[1] + drjl[0]*fj[1] + drkl[0]*fk[1];
  v[4] = dril[0]*fi[2] + drjl[0]*fj[2] + drkl[0]*fk[2];
  v[5] = dril[1]*fi[2] + drjl[1]*fj[2] + drkl[1]*fk[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    a_vatom(i,0) += 0.25 * v[0]; a_vatom(i,1) += 0.25 * v[1]; a_vatom(i,2) += 0.25 * v[2];
    a_vatom(i,3) += 0.25 * v[3]; a_vatom(i,4) += 0.25 * v[4]; a_vatom(i,5) += 0.25 * v[5];
    a_vatom(j,0) += 0.25 * v[0]; a_vatom(j,1) += 0.25 * v[1]; a_vatom(j,2) += 0.25 * v[2];
    a_vatom(j,3) += 0.25 * v[3]; a_vatom(j,4) += 0.25 * v[4]; a_vatom(j,5) += 0.25 * v[5];
    a_vatom(k,0) += 0.25 * v[0]; a_vatom(k,1) += 0.25 * v[1]; a_vatom(k,2) += 0.25 * v[2];
    a_vatom(k,3) += 0.25 * v[3]; a_vatom(k,4) += 0.25 * v[4]; a_vatom(k,5) += 0.25 * v[5];
    a_vatom(l,0) += 0.25 * v[0]; a_vatom(l,1) += 0.25 * v[1]; a_vatom(l,2) += 0.25 * v[2];
    a_vatom(l,3) += 0.25 * v[3]; a_vatom(l,4) += 0.25 * v[4]; a_vatom(l,5) += 0.25 * v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally3_atom(EV_FLOAT_REAX &ev, const int &i, const int & /*j*/,
                                                const int & /*k*/, F_FLOAT *fj, F_FLOAT *fk,
                                                F_FLOAT *drji, F_FLOAT *drjk) const
{
  F_FLOAT v[6];

  v[0] = THIRD * (drji[0]*fj[0] + drjk[0]*fk[0]);
  v[1] = THIRD * (drji[1]*fj[1] + drjk[1]*fk[1]);
  v[2] = THIRD * (drji[2]*fj[2] + drjk[2]*fk[2]);
  v[3] = THIRD * (drji[0]*fj[1] + drjk[0]*fk[1]);
  v[4] = THIRD * (drji[0]*fj[2] + drjk[0]*fk[2]);
  v[5] = THIRD * (drji[1]*fj[2] + drjk[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    d_vatom(i,0) += v[0]; d_vatom(i,1) += v[1]; d_vatom(i,2) += v[2];
    d_vatom(i,3) += v[3]; d_vatom(i,4) += v[4]; d_vatom(i,5) += v[5];
  }
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag and vflag
   see pair::ev_setup() for values of eflag_* and vflag_*
   VIRIAL_CENTROID bitflag is not yet supported by ReaxFF
------------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::ev_setup(int eflag, int vflag, int)
{
  int i;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_either = vflag;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  vflag_atom = vflag & VIRIAL_ATOM;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  // zero accumulators

  if (eflag_global) eng_vdwl = eng_coul = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) Kokkos::deep_copy(d_eatom,0.0);
  if (vflag_atom) Kokkos::deep_copy(d_vatom,0.0);

  // if vflag_global = VIRIAL_FDOTR and pair::compute() calls virial_fdotr_compute()
  // compute global virial via (F dot r) instead of via pairwise summation
  // unset other flags as appropriate

  if (vflag_global == VIRIAL_FDOTR && no_virial_fdotr_compute == 0) {
    vflag_fdotr = 1;
    vflag_global = 0;
    if (vflag_atom == 0) vflag_either = 0;
    if (vflag_either == 0 && eflag_either == 0) evflag = 0;
  } else vflag_fdotr = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = 0.0;

  if (cut_hbsq > 0.0) {
    bytes += (double)nmax*3*sizeof(int);
    bytes += (double)maxhb*nmax*sizeof(int);
  }
  bytes += (double)nmax*2*sizeof(int);
  bytes += (double)maxbo*nmax*sizeof(int);

  bytes += (double)nmax*17*sizeof(F_FLOAT);
  bytes += (double)maxbo*nmax*34*sizeof(F_FLOAT);

  // FixReaxFFSpecies
  if (fixspecies_flag) {
    bytes += (double)MAXSPECBOND*nmax*sizeof(tagint);
    bytes += (double)MAXSPECBOND*nmax*sizeof(F_FLOAT);
  }

  // FixReaxFFBonds
  bytes += (double)maxbo*nmax*sizeof(tagint);
  bytes += (double)maxbo*nmax*sizeof(F_FLOAT);
  bytes += (double)nmax*sizeof(int);

  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::FindBond(int &numbonds)
{
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondZero>(0,nmax),*this);

  bo_cut_bond = api->control->bg_cut;

  atomKK->sync(execution_space,TAG_MASK);
  tag = atomKK->k_tag.view<DeviceType>();

  const int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;

  numbonds = 0;
  PairReaxKokkosFindBondFunctor<DeviceType> find_bond_functor(this);
  Kokkos::parallel_reduce(inum,find_bond_functor,numbonds);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondZero, const int &i) const {
  d_numneigh_bonds[i] = 0;
  for (int j = 0; j < maxbo; j++) {
    d_neighid(i,j) = 0;
    d_abo(i,j) = 0.0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::calculate_find_bond_item(int ii, int &numbonds) const
{
  const int i = d_ilist[ii];
  int nj = 0;

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];
  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    const tagint jtag = tag[j];
    const int j_index = jj - j_start;
    double bo_tmp = d_BO(i,j_index);

    if (bo_tmp > bo_cut_bond) {
      d_neighid(i,nj) = jtag;
      d_abo(i,nj) = bo_tmp;
      nj++;
    }
  }
  d_numneigh_bonds[i] = nj;
  if (nj > numbonds) numbonds = nj;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::PackBondBuffer(DAT::tdual_ffloat_1d k_buf, int &nbuf_local)
{
  d_buf = k_buf.view<DeviceType>();
  k_params_sing.template sync<DeviceType>();
  atomKK->sync(execution_space,TAG_MASK|TYPE_MASK|Q_MASK|MOLECULE_MASK);

  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  if (atom->molecule)
    molecule = atomKK->k_molecule.view<DeviceType>();

  copymode = 1;
  nlocal = atomKK->nlocal;
  PairReaxKokkosPackBondBufferFunctor<DeviceType> pack_bond_buffer_functor(this);
  Kokkos::parallel_scan(nlocal,pack_bond_buffer_functor);
  copymode = 0;

  k_buf.modify<DeviceType>();
  k_nbuf_local.modify<DeviceType>();

  k_buf.sync<LMPHostType>();
  k_nbuf_local.sync<LMPHostType>();
  nbuf_local = k_nbuf_local.h_view();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::pack_bond_buffer_item(int i, int &j, const bool &final) const
{
  if (i == 0)
    j += 2;

  if (final) {
    d_buf[j-1] = tag[i];
    d_buf[j+0] = type[i];
    d_buf[j+1] = d_total_bo[i];
    d_buf[j+2] = paramssing(type[i]).nlp_opt - d_Delta_lp[i];
    d_buf[j+3] = q[i];
    d_buf[j+4] = d_numneigh_bonds[i];
  }
  const int numbonds = d_numneigh_bonds[i];

  if (final) {
    for (int k = 5; k < 5+numbonds; k++) {
      d_buf[j+k] = d_neighid(i,k-5);
    }
  }
  j += (5+numbonds);

  if (final) {
    if (!molecule.data()) d_buf[j] = 0.0;
    else d_buf[j] = molecule[i];
  }
  j++;

  if (final) {
    for (int k = 0; k < numbonds; k++) {
      d_buf[j+k] = d_abo(i,k);
    }
  }
  j += (1+numbonds);

  if (final && i == nlocal-1)
    k_nbuf_local.view<DeviceType>()() = j - 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::FindBondSpecies()
{

  if (nmax > (int)k_tmpid.extent(0)) {
    memoryKK->destroy_kokkos(k_tmpid,tmpid);
    memoryKK->destroy_kokkos(k_tmpbo,tmpbo);
    memoryKK->create_kokkos(k_tmpid,tmpid,nmax,MAXSPECBOND,"pair:tmpid");
    memoryKK->create_kokkos(k_tmpbo,tmpbo,nmax,MAXSPECBOND,"pair:tmpbo");
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondSpeciesZero>(0,nmax),*this);

  nlocal = atomKK->nlocal;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondSpecies>(0,nlocal),*this);
  copymode = 0;

  // NOTE: Could improve performance if a Kokkos version of ComputeSpecAtom is added

  k_tmpbo.modify<DeviceType>();
  k_tmpid.modify<DeviceType>();
  k_error_flag.modify<DeviceType>();

  k_tmpbo.sync<LMPHostType>();
  k_tmpid.sync<LMPHostType>();
  k_error_flag.sync<LMPHostType>();

  if (k_error_flag.h_view())
    error->all(FLERR,"Increase MAXSPECBOND in reaxff_defs.h");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondSpeciesZero, const int &i) const {
  for (int j = 0; j < MAXSPECBOND; j++) {
    k_tmpbo.view<DeviceType>()(i,j) = 0.0;
    k_tmpid.view<DeviceType>()(i,j) = 0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondSpecies, const int &i) const {
  int nj = 0;

  const int j_start = d_bo_first[i];
  const int j_end = j_start + d_bo_num[i];
  for (int jj = j_start; jj < j_end; jj++) {
    int j = d_bo_list[jj];
    j &= NEIGHMASK;
    if (j < i) continue;
    const int j_index = jj - j_start;

    double bo_tmp = d_BO(i,j_index);

    if (bo_tmp >= 0.10) { // Why is this a hardcoded value?
      k_tmpid.view<DeviceType>()(i,nj) = j;
      k_tmpbo.view<DeviceType>()(i,nj) = bo_tmp;
      nj++;
      if (nj > MAXSPECBOND) k_error_flag.view<DeviceType>()() = 1;
    }
  }
}

template class PairReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairReaxFFKokkos<LMPHostType>;
#endif
}
