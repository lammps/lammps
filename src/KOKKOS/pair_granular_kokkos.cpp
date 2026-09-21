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
   Kokkos version of the modular pair style granular.

   The contact physics is not reimplemented here: the kernel calls the same
   stateless functions in src/GRANULAR/gran_sub_mod_*_kernel.h that the host
   sub-model classes call, so the two code paths cannot drift apart.  Only
   coefficient parsing, mixing and restart handling stay on the host.

   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_granular_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_neigh_history_kokkos.h"
#include "force.h"
#include "gran_sub_mod.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_heat.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_rolling.h"
#include "gran_sub_mod_tangential.h"
#include "gran_sub_mod_twisting.h"
#include "granular_model.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace LAMMPS_NS::Granular_NS;
using namespace LAMMPS_NS::Granular_NS::GranKernel;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranularKokkos<DeviceType>::PairGranularKokkos(LAMMPS *lmp) : PairGranular(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK |
    TYPE_MASK | MASK_MASK | ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  fix_historyKK = nullptr;
  use_history_kk = 0;
  size_history_kk = 0;
  extra_models = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranularKokkos<DeviceType>::~PairGranularKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ----------------------------------------------------------------------
   ask the base class for the Kokkos flavor of the neighbor history fix
------------------------------------------------------------------------- */

template<class DeviceType>
std::string PairGranularKokkos<DeviceType>::history_fix_command()
{
  std::string style = (execution_space == Device) ? "NEIGH_HISTORY/KK/DEVICE"
                                                  : "NEIGH_HISTORY/KK/HOST";
  return fmt::format("{} all {} {}", id_history, style, size_history);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGranularKokkos<DeviceType>::init_style()
{
  // PairGranular::init_style() creates the history fix through the
  // history_fix_command() hook overridden above

  PairGranular::init_style();

  if (use_history) {
    fix_historyKK = dynamic_cast<FixNeighHistoryKokkos<DeviceType> *>(fix_history);
    if (!fix_historyKK)
      error->all(FLERR,"Pair granular/kk could not create the Kokkos neighbor history fix");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Pair granular/kk requires a half or half/thread neighbor list");

  // fix neigh/history/kk only implements the newton off exchange path

  if (use_history && force->newton_pair)
    error->all(FLERR,"Pair granular/kk with contact history requires newton off");

  if (fix_rigid)
    error->all(FLERR,"Pair granular/kk does not yet support rigid body masses");

  if (size_history > GRAN_KK_MAX_HISTORY)
    error->all(FLERR,"Pair granular/kk does not support a contact history of {} values",
               size_history);

  use_history_kk = use_history;
  size_history_kk = size_history;
}

/* ----------------------------------------------------------------------
   translate the host sub-model objects into plain data for the kernel
   run once per setup, after coefficient mixing has produced all models
------------------------------------------------------------------------- */

template<class DeviceType>
void PairGranularKokkos<DeviceType>::pack_models()
{
  auto h_models = Kokkos::View<GranularModelKK<KK_FLOAT>*,Kokkos::HostSpace>(
      "pair_granular:models_host",nmodels);
  auto h_type_model = Kokkos::View<int**,Kokkos::LayoutRight,Kokkos::HostSpace>(
      "pair_granular:type_model_host",atom->ntypes+1,atom->ntypes+1);

  extra_models = 0;

  for (int n = 0; n < nmodels; n++) {
    auto *gm = models_list[n];
    auto &m = h_models(n);
    memset(&m,0,sizeof(GranularModelKK<KK_FLOAT>));

    // ---- normal ----

    const std::string &normal = gm->normal_model->name;
    if (normal == "hooke") m.normal.model = GRAN_NORMAL_HOOKE;
    else if (normal == "hertz") m.normal.model = GRAN_NORMAL_HERTZ;
    else if (normal == "hertz/material") m.normal.model = GRAN_NORMAL_HERTZ_MATERIAL;
    else if (normal == "dmt") m.normal.model = GRAN_NORMAL_DMT;
    else if (normal == "jkr") m.normal.model = GRAN_NORMAL_JKR;
    else error->all(FLERR,"Pair granular/kk does not support normal model {}",normal);

    GranNormalParams<double> pn;
    gm->normal_model->fill_kernel_params(pn);
    m.normal.k = pn.k;
    m.normal.cohesion = pn.cohesion;
    m.normal.Emix = pn.Emix;

    // ---- damping ----

    const std::string &damping = gm->damping_model->name;
    if (damping == "none") m.damping.model = GRAN_DAMPING_NONE;
    else if (damping == "velocity") m.damping.model = GRAN_DAMPING_VELOCITY;
    else if (damping == "mass_velocity") m.damping.model = GRAN_DAMPING_MASS_VELOCITY;
    else if (damping == "viscoelastic") m.damping.model = GRAN_DAMPING_VISCOELASTIC;
    else if ((damping == "tsuji") || (damping == "coeff_restitution"))
      m.damping.model = GRAN_DAMPING_TSUJI;
    else error->all(FLERR,"Pair granular/kk does not support damping model {}",damping);

    GranDampingParams<double> pd;
    gm->damping_model->fill_kernel_params(pd);
    m.damping.damp = pd.damp;

    // ---- tangential ----

    const std::string &tangential = gm->tangential_model->name;
    if (tangential == "linear_nohistory")
      m.tangential.model = GRAN_TANGENTIAL_LINEAR_NOHISTORY;
    else if (tangential == "linear_history")
      m.tangential.model = GRAN_TANGENTIAL_LINEAR_HISTORY;
    else if (tangential == "linear_history_classic")
      m.tangential.model = GRAN_TANGENTIAL_LINEAR_HISTORY_CLASSIC;
    else if (tangential == "mindlin_classic")
      m.tangential.model = GRAN_TANGENTIAL_MINDLIN_CLASSIC;
    else if (utils::strmatch(tangential,"^mindlin"))
      m.tangential.model = GRAN_TANGENTIAL_MINDLIN;
    else error->all(FLERR,"Pair granular/kk does not support tangential model {}",tangential);

    GranTangentialParams<double> pt;
    gm->tangential_model->fill_kernel_params(pt);
    m.tangential.k = pt.k;
    m.tangential.xt = pt.xt;
    m.tangential.mu = pt.mu;
    m.tangential.mindlin_force = pt.mindlin_force;
    m.tangential.mindlin_rescale = pt.mindlin_rescale;
    m.tangential.contact_radius_flag = pt.contact_radius_flag;

    // ---- rolling and twisting ----

    const std::string &rolling = gm->rolling_model->name;
    if (rolling == "none") m.rolling.model = GRAN_ROLLING_NONE;
    else if (rolling == "sds") m.rolling.model = GRAN_ROLLING_SDS;
    else error->all(FLERR,"Pair granular/kk does not support rolling model {}",rolling);
    GranRollingParams<double> pr;
    gm->rolling_model->fill_kernel_params(pr);
    m.rolling.k = pr.k;
    m.rolling.gamma = pr.gamma;
    m.rolling.mu = pr.mu;

    const std::string &twisting = gm->twisting_model->name;
    if (twisting == "none") m.twisting.model = GRAN_TWISTING_NONE;
    else if (twisting == "marshall") m.twisting.model = GRAN_TWISTING_MARSHALL;
    else if (twisting == "sds") m.twisting.model = GRAN_TWISTING_SDS;
    else error->all(FLERR,"Pair granular/kk does not support twisting model {}",twisting);
    GranTwistingParams<double> pw;
    gm->twisting_model->fill_kernel_params(pw);
    m.twisting.k = pw.k;
    m.twisting.damp = pw.damp;
    m.twisting.mu = pw.mu;
    m.twisting.k_tang = pw.k_tang;
    m.twisting.mu_tang = pw.mu_tang;

    if ((m.rolling.model != GRAN_ROLLING_NONE) || (m.twisting.model != GRAN_TWISTING_NONE))
      extra_models = 1;

    // ---- unsupported options ----

    if (gm->heat_model->name != "none")
      error->all(FLERR,"Pair granular/kk does not yet support heat conduction models");
    if (gm->synchronized_verlet)
      error->all(FLERR,"Pair granular/kk does not yet support synchronized verlet");
    if (gm->nondefault_history_transfer)
      error->all(FLERR,"Pair granular/kk does not yet support the {} tangential model "
                 "because fix neigh/history/kk cannot transfer its history",tangential);

    // ---- history layout ----

    m.tangential_index = gm->tangential_model->history_index;
    m.rolling_index = gm->rolling_model->history_index;
    m.twisting_index = gm->twisting_model->history_index;
    if (m.tangential_index != 0)
      error->all(FLERR,"Pair granular/kk expects the tangential history to come first");

    m.tangential_size = gm->tangential_model->size_history;
    m.rolling_size = gm->rolling_model->size_history;
    m.twisting_size = gm->twisting_model->size_history;
    if ((m.tangential_size > GRAN_KK_TANGENTIAL_HISTORY) ||
        (m.rolling_size > GRAN_KK_ROLLING_HISTORY) ||
        (m.twisting_size > GRAN_KK_TWISTING_HISTORY))
      error->all(FLERR,"Pair granular/kk does not support this contact history size");

    m.limit_damping = gm->limit_damping;
    m.contact_radius_flag = gm->get_contact_radius_flag();
  }

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++)
      h_type_model(i,j) = types_indices[i][j];

  d_models = Kokkos::View<GranularModelKK<KK_FLOAT>*,DeviceType>(
      "pair_granular:models",nmodels);
  d_type_model = Kokkos::View<int**,Kokkos::LayoutRight,DeviceType>(
      "pair_granular:type_model",atom->ntypes+1,atom->ntypes+1);
  Kokkos::deep_copy(d_models,h_models);
  Kokkos::deep_copy(d_type_model,h_type_model);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGranularKokkos<DeviceType>::setup()
{
  PairGranular::setup();
  pack_models();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairGranularKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // Verlet::setup_minimal() (run ... pre no) skips Force::setup(), so make
  // sure the model tables exist before the first kernel launch

  if (d_models.extent(0) != (size_t) nmodels) pack_models();

  history_update = (update->setupflag == 0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  dt_kk = static_cast<KK_FLOAT>(update->dt);
  for (int k = 0; k < 4; k++) special_lj[k] = static_cast<KK_FLOAT>(force->special_lj[k]);

  const int inum = list->inum;
  auto *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  if (use_history_kk) {
    fix_historyKK->k_firstflag.template sync<DeviceType>();
    fix_historyKK->k_firstvalue.template sync<DeviceType>();
    d_firsttouch = fix_historyKK->k_firstflag.template view<DeviceType>();
    d_firsthistory = fix_historyKK->k_firstvalue.template view<DeviceType>();
  }

  EV_FLOAT ev;

#define LMP_GRAN_KK_RUN(NF,NP,EX)                                                             \
  if (vflag_either)                                                                           \
    Kokkos::parallel_reduce(                                                                  \
      Kokkos::RangePolicy<DeviceType,TagPairGranularCompute<NF,NP,1,EX>>(0,inum),*this,ev);    \
  else                                                                                        \
    Kokkos::parallel_for(                                                                     \
      Kokkos::RangePolicy<DeviceType,TagPairGranularCompute<NF,NP,0,EX>>(0,inum),*this)

#define LMP_GRAN_KK_RUN_EXTRA(NF,NP)                                                          \
  if (extra_models) { LMP_GRAN_KK_RUN(NF,NP,1); } else { LMP_GRAN_KK_RUN(NF,NP,0); }

  if (neighflag == HALF) {
    if (newton_pair) { LMP_GRAN_KK_RUN_EXTRA(HALF,1); }
    else             { LMP_GRAN_KK_RUN_EXTRA(HALF,0); }
  } else {
    if (newton_pair) { LMP_GRAN_KK_RUN_EXTRA(HALFTHREAD,1); }
    else             { LMP_GRAN_KK_RUN_EXTRA(HALFTHREAD,0); }
  }

#undef LMP_GRAN_KK_RUN_EXTRA
#undef LMP_GRAN_KK_RUN

  if (use_history_kk) {
    fix_historyKK->k_firstflag.template modify<DeviceType>();
    fix_historyKK->k_firstvalue.template modify<DeviceType>();
  }

  if (vflag_global) for (int k = 0; k < 6; k++) virial[k] += (double) ev.v[k];
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   one contact; mirrors GranularModel::check_contact() and
   GranularModel::calculate_forces()
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int EXTRA>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::operator()(
    TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,EXTRA>, const int ii, EV_FLOAT &ev) const
{
  using AtomicView = Kokkos::View<KK_ACC_FLOAT*[3],typename DAT::t_kkacc_1d_3::array_layout,
    typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>>;
  AtomicView a_f = f;
  AtomicView a_torque = torque;

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0), ytmp = x(i,1), ztmp = x(i,2);
  const KK_FLOAT vxi = v(i,0), vyi = v(i,1), vzi = v(i,2);
  const KK_FLOAT wxi = omega(i,0), wyi = omega(i,1), wzi = omega(i,2);
  const KK_FLOAT radi = radius[i];
  const KK_FLOAT mi_mass = rmass[i];
  const int itype = type[i];
  const int maski = mask[i];
  const int jnum = d_numneigh[i];
  const int nhist = size_history_kk;

  KK_ACC_FLOAT fxi = 0.0, fyi = 0.0, fzi = 0.0;
  KK_ACC_FLOAT txi = 0.0, tyi = 0.0, tzi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;
    if (factor_lj == 0.0) continue;

    const GranularModelKK<KK_FLOAT> &m = d_models(d_type_model(itype,type[j]));

    // ---- check_contact() ----

    const KK_FLOAT dx[3] = {xtmp - x(j,0), ytmp - x(j,1), ztmp - x(j,2)};
    const KK_FLOAT rsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const KK_FLOAT radj = radius[j];
    const KK_FLOAT radsum = radi + radj;
    const KK_FLOAT Reff = radi * radj / radsum;

    // only JKR reads the stored contact flag; skip the load for the others
    const bool touch_prev = use_history_kk && (m.normal.model == GRAN_NORMAL_JKR) &&
      (d_firsttouch(i,jj) != 0);
    const bool touchflag = gran_normal_touch(m.normal,rsq,radsum,Reff,touch_prev);

    if (!touchflag) {
      if (use_history_kk) {
        d_firsttouch(i,jj) = 0;
        for (int k = 0; k < nhist; k++) d_firsthistory(i,nhist*jj+k) = 0.0;
      }
      continue;
    }
    if (use_history_kk) d_firsttouch(i,jj) = 1;

    // ---- effective mass ----

    const KK_FLOAT mj_mass = rmass[j];
    KK_FLOAT meff = mi_mass * mj_mass / (mi_mass + mj_mass);
    if (maski & freeze_group_bit) meff = mj_mass;
    if (mask[j] & freeze_group_bit) meff = mi_mass;

    // ---- geometry, GranularModel::calculate_forces() ----

    const KK_FLOAT r = GRAN_MATH::sqrt(rsq);
    const KK_FLOAT rinv = 1.0 / r;
    const KK_FLOAT delta = radsum - r;
    const KK_FLOAT dR = delta * Reff;

    const KK_FLOAT vr[3] = {vxi - v(j,0), vyi - v(j,1), vzi - v(j,2)};
    const KK_FLOAT nx[3] = {rinv*dx[0], rinv*dx[1], rinv*dx[2]};

    const KK_FLOAT vnnr = vr[0]*nx[0] + vr[1]*nx[1] + vr[2]*nx[2];
    const KK_FLOAT vt[3] = {vr[0] - vnnr*nx[0], vr[1] - vnnr*nx[1], vr[2] - vnnr*nx[2]};

    const KK_FLOAT wr[3] = {radi*wxi + radj*omega(j,0),
                            radi*wyi + radj*omega(j,1),
                            radi*wzi + radj*omega(j,2)};
    KK_FLOAT tmp3[3];
    gk_cross3(wr,nx,tmp3);
    const KK_FLOAT vtr[3] = {vt[0]-tmp3[0], vt[1]-tmp3[1], vt[2]-tmp3[2]};
    const KK_FLOAT vrel = gk_len3(vtr);

    // ---- normal and damping ----

    KK_FLOAT contact_radius = 0.0;
    if (m.contact_radius_flag)
      contact_radius = gran_normal_contact_radius(m.normal,Reff,dR);

    KK_FLOAT Fne, F_pulloff;
    const KK_FLOAT Fnormal = gran_normal_force(m.normal,contact_radius,delta,Reff,Fne,F_pulloff);
    const KK_FLOAT damp_prefactor =
      gran_damping_prefactor(m.damping,meff,contact_radius,Fnormal,delta);

    KK_FLOAT Fntot = Fnormal - damp_prefactor * vnnr;
    if (m.limit_damping && (Fntot < 0.0)) Fntot = 0.0;

    const KK_FLOAT Fncrit = gran_normal_fncrit(m.normal,Fntot,Fne,F_pulloff);

    // ---- tangential ----

    GranTangentialState<KK_FLOAT> ts;
    ts.nx = nx;
    ts.nx_unrotated = nx;
    ts.vtr = vtr;
    ts.vrel = vrel;
    ts.dt = dt_kk;
    ts.contact_radius = contact_radius;
    ts.damp_prefactor = damp_prefactor;
    ts.Fncrit = Fncrit;
    ts.synchronized_verlet = 0;
    ts.history_update = history_update;

    KK_FLOAT fs[3] = {0.0, 0.0, 0.0};

    // one fixed size array per sub-model, indexed by constants only, so
    // the history stays in registers; tangential_index is 0 for every
    // supported model, checked in pack_models()

    const int hbase = nhist * jj;
    KK_FLOAT ht[GRAN_KK_TANGENTIAL_HISTORY] = {0.0};
    if (use_history_kk) {
      for (int k = 0; k < GRAN_KK_TANGENTIAL_HISTORY; k++)
        if (k < m.tangential_size) ht[k] = d_firsthistory(i,hbase+k);
    }

    gran_tangential_forces(m.tangential,ts,ht,fs);

    // the Coulomb rescale writes history even when history_update is 0
    if (use_history_kk) {
      for (int k = 0; k < GRAN_KK_TANGENTIAL_HISTORY; k++)
        if (k < m.tangential_size) d_firsthistory(i,hbase+k) = ht[k];
    }

    // ---- sum normal + tangential ----

    KK_FLOAT forces[3] = {Fntot*nx[0] + fs[0], Fntot*nx[1] + fs[1], Fntot*nx[2] + fs[2]};

    KK_FLOAT torquesi[3], torquesj[3];
    gk_cross3(nx,fs,torquesi);
    gk_scale3((KK_FLOAT) -1.0,torquesi);
    gk_copy3(torquesi,torquesj);
    gk_scale3(radi - (KK_FLOAT) 0.5*delta,torquesi);
    gk_scale3(radj - (KK_FLOAT) 0.5*delta,torquesj);

    // ---- rolling and twisting ----

    if (EXTRA) {
      const KK_FLOAT relrot[3] = {wxi - omega(j,0), wyi - omega(j,1), wzi - omega(j,2)};

      if (m.rolling.model != GRAN_ROLLING_NONE) {
        const KK_FLOAT vrl[3] = {Reff * (relrot[1]*nx[2] - relrot[2]*nx[1]),
                                 Reff * (relrot[2]*nx[0] - relrot[0]*nx[2]),
                                 Reff * (relrot[0]*nx[1] - relrot[1]*nx[0])};
        GranRollingState<KK_FLOAT> rs;
        rs.nx = nx;
        rs.nx_unrotated = nx;
        rs.vrl = vrl;
        rs.dt = dt_kk;
        rs.Fncrit = Fncrit;
        rs.synchronized_verlet = 0;
        rs.history_update = history_update;

        KK_FLOAT hr[GRAN_KK_ROLLING_HISTORY] = {0.0};
        const int hroll = hbase + m.rolling_index;
        for (int k = 0; k < GRAN_KK_ROLLING_HISTORY; k++)
          if (k < m.rolling_size) hr[k] = d_firsthistory(i,hroll+k);

        KK_FLOAT fr[3];
        gran_rolling_forces(m.rolling,rs,hr,fr);

        for (int k = 0; k < GRAN_KK_ROLLING_HISTORY; k++)
          if (k < m.rolling_size) d_firsthistory(i,hroll+k) = hr[k];

        KK_FLOAT torroll[3];
        gk_cross3(nx,fr,torroll);
        gk_scale3(Reff,torroll);
        gk_add3(torquesi,torroll,torquesi);
        gk_sub3(torquesj,torroll,torquesj);
      }

      if (m.twisting.model != GRAN_TWISTING_NONE) {
        const KK_FLOAT magtwist = relrot[0]*nx[0] + relrot[1]*nx[1] + relrot[2]*nx[2];
        const KK_FLOAT tangential_damp = gran_tangential_damp(m.tangential,damp_prefactor);

        KK_FLOAT hw[GRAN_KK_TWISTING_HISTORY] = {0.0};
        const int htwist = hbase + m.twisting_index;
        for (int k = 0; k < GRAN_KK_TWISTING_HISTORY; k++)
          if (k < m.twisting_size) hw[k] = d_firsthistory(i,htwist+k);

        const KK_FLOAT magtortwist =
          gran_twisting_forces(m.twisting,magtwist,dt_kk,contact_radius,Fncrit,
                               tangential_damp,history_update,hw);

        for (int k = 0; k < GRAN_KK_TWISTING_HISTORY; k++)
          if (k < m.twisting_size) d_firsthistory(i,htwist+k) = hw[k];

        KK_FLOAT tortwist[3];
        gk_scale3(magtortwist,nx,tortwist);
        gk_add3(torquesi,tortwist,torquesi);
        gk_sub3(torquesj,tortwist,torquesj);
      }
    }

    // ---- apply ----

    gk_scale3(factor_lj,forces);
    fxi += forces[0]; fyi += forces[1]; fzi += forces[2];

    gk_scale3(factor_lj,torquesi);
    txi += torquesi[0]; tyi += torquesi[1]; tzi += torquesi[2];

    if (NEWTON_PAIR || (j < nlocal)) {
      a_f(j,0) -= forces[0];
      a_f(j,1) -= forces[1];
      a_f(j,2) -= forces[2];
      gk_scale3(factor_lj,torquesj);
      a_torque(j,0) += torquesj[0];
      a_torque(j,1) += torquesj[1];
      a_torque(j,2) += torquesj[2];
    }

    if (VFLAG)
      ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,forces[0],forces[1],forces[2],
                                          dx[0],dx[1],dx[2]);
  }

  a_f(i,0) += fxi; a_f(i,1) += fyi; a_f(i,2) += fzi;
  a_torque(i,0) += txi; a_torque(i,1) += tyi; a_torque(i,2) += tzi;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int EXTRA>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::operator()(
    TagPairGranularCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,EXTRA> tag, const int ii) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,VFLAG,EXTRA>(tag,ii,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz, KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  Kokkos::View<KK_ACC_FLOAT*[6],typename DAT::t_kkacc_1d_6::array_layout,
    typename KKDevice<DeviceType>::value,
    Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> v_vatom = d_vatom;

  if (vflag_global || vflag_atom) {
    const KK_ACC_FLOAT v_acc[6] =
      { static_cast<KK_ACC_FLOAT>(delx*fx),
        static_cast<KK_ACC_FLOAT>(dely*fy),
        static_cast<KK_ACC_FLOAT>(delz*fz),
        static_cast<KK_ACC_FLOAT>(delx*fy),
        static_cast<KK_ACC_FLOAT>(delx*fz),
        static_cast<KK_ACC_FLOAT>(dely*fz) };

    if (vflag_global) {
      if (NEWTON_PAIR) {
        for (int n = 0; n < 6; n++) ev.v[n] += v_acc[n];
      } else {
        if (i < nlocal)
          for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
        if (j < nlocal)
          for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
      }
    }

    if (vflag_atom) {
      if (NEWTON_PAIR || (i < nlocal))
        for (int n = 0; n < 6; n++) v_vatom(i,n) += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
      if (NEWTON_PAIR || (j < nlocal))
        for (int n = 0; n < 6; n++) v_vatom(j,n) += static_cast<KK_ACC_FLOAT>(0.5) * v_acc[n];
    }
  }
}

namespace LAMMPS_NS {
template class PairGranularKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGranularKokkos<LMPHostType>;
#endif
}
