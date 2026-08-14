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
   Contributing Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_switch3_coulgauss_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJSwitch3CoulGaussLongKokkos<DeviceType>::PairLJSwitch3CoulGaussLongKokkos(LAMMPS *lmp)
  : PairLJSwitch3CoulGaussLong(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJSwitch3CoulGaussLongKokkos<DeviceType>::~PairLJSwitch3CoulGaussLongKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_ljsq,cut_ljsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSwitch3CoulGaussLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

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
  k_cutsq.template sync<DeviceType>();
  k_cut_ljsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;
  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  truncw_kk = static_cast<KK_FLOAT>(truncw);
  truncwi_kk = static_cast<KK_FLOAT>(truncwi);

  EV_FLOAT ev;

  copymode = 1;

  ev = pair_compute<PairLJSwitch3CoulGaussLongKokkos<DeviceType>,void>
    (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   LJ 12-6 pair force with Switch3 truncation
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSwitch3CoulGaussLongKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT lj3_val = STACKPARAMS ? m_params[itype][jtype].lj3 : params(itype,jtype).lj3;
  const KK_FLOAT lj4_val = STACKPARAMS ? m_params[itype][jtype].lj4 : params(itype,jtype).lj4;
  const KK_FLOAT offset_val = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;
  const KK_FLOAT cut_lj_val = STACKPARAMS ? m_params[itype][jtype].cut_lj : params(itype,jtype).cut_lj;

  KK_FLOAT forcelj = r6inv*(static_cast<KK_FLOAT>(12.0)*lj3_val*r6inv
                            - static_cast<KK_FLOAT>(6.0)*lj4_val);
  KK_FLOAT evdwl = r6inv*(lj3_val*r6inv - lj4_val) - offset_val;

  // Switch3 truncation: smooth to zero near cutoff
  if (truncw_kk > static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT r = Kokkos::sqrt(rsq);
    if (r > cut_lj_val - truncw_kk) {
      const KK_FLOAT trx = (cut_lj_val - r) * truncwi_kk;
      const KK_FLOAT tr = trx*trx*(static_cast<KK_FLOAT>(3.0) - static_cast<KK_FLOAT>(2.0)*trx);
      const KK_FLOAT ftr = static_cast<KK_FLOAT>(6.0)*trx*(static_cast<KK_FLOAT>(1.0)-trx)*r*truncwi_kk;
      forcelj = forcelj*tr + evdwl*ftr;
    }
  }

  return forcelj * r2inv;
}

/* ----------------------------------------------------------------------
   Ewald long-range Coulomb force + Gaussian charge correction
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSwitch3CoulGaussLongKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc1 = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                         t * (static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                         t * static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) * rinv;

  // Ewald Coulomb force
  KK_FLOAT forcecoul = prefactor * (erfc1 + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  // Gaussian correction (only within VdW cutoff)
  const KK_FLOAT cut_ljsq_val = STACKPARAMS ? m_params[itype][jtype].cut_ljsq : params(itype,jtype).cut_ljsq;
  const KK_FLOAT lj2_val = STACKPARAMS ? m_params[itype][jtype].lj2 : params(itype,jtype).lj2;
  if (rsq < cut_ljsq_val && lj2_val > static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT rrij = lj2_val * r;
    const KK_FLOAT expn2 = Kokkos::exp(-rrij*rrij);
    const KK_FLOAT erfc2 = Kokkos::erfc(rrij);
    const KK_FLOAT prefactor2 = -qqrd2e * qtmp * q(j) * rinv;
    forcecoul += factor_coul * prefactor2 * (erfc2 + static_cast<KK_FLOAT>(EWALD_F)*rrij*expn2);
  }

  return forcecoul * r2inv;
}

/* ----------------------------------------------------------------------
   LJ 12-6 pair potential energy with Switch3 truncation
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSwitch3CoulGaussLongKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
               const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT lj3_val = STACKPARAMS ? m_params[itype][jtype].lj3 : params(itype,jtype).lj3;
  const KK_FLOAT lj4_val = STACKPARAMS ? m_params[itype][jtype].lj4 : params(itype,jtype).lj4;
  const KK_FLOAT offset_val = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;
  const KK_FLOAT cut_lj_val = STACKPARAMS ? m_params[itype][jtype].cut_lj : params(itype,jtype).cut_lj;

  KK_FLOAT evdwl = r6inv*(lj3_val*r6inv - lj4_val) - offset_val;

  // Switch3 truncation
  if (truncw_kk > static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT r = Kokkos::sqrt(rsq);
    if (r > cut_lj_val - truncw_kk) {
      const KK_FLOAT trx = (cut_lj_val - r) * truncwi_kk;
      const KK_FLOAT tr = trx*trx*(static_cast<KK_FLOAT>(3.0) - static_cast<KK_FLOAT>(2.0)*trx);
      evdwl *= tr;
    }
  }

  return evdwl;
}

/* ----------------------------------------------------------------------
   Ewald long-range Coulomb energy + Gaussian charge correction
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSwitch3CoulGaussLongKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
               const int& itype, const int& jtype,
               const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc1 = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                         t * (static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                         t * static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) * rinv;

  KK_FLOAT ecoul = prefactor * erfc1;
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  // Gaussian correction (only within VdW cutoff)
  const KK_FLOAT cut_ljsq_val = STACKPARAMS ? m_params[itype][jtype].cut_ljsq : params(itype,jtype).cut_ljsq;
  const KK_FLOAT lj2_val = STACKPARAMS ? m_params[itype][jtype].lj2 : params(itype,jtype).lj2;
  if (rsq < cut_ljsq_val && lj2_val > static_cast<KK_FLOAT>(0.0)) {
    const KK_FLOAT rrij = lj2_val * r;
    const KK_FLOAT erfc2 = Kokkos::erfc(rrij);
    const KK_FLOAT prefactor2 = -qqrd2e * qtmp * q(j) * rinv;
    ecoul += factor_coul * prefactor2 * erfc2;
  }

  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSwitch3CoulGaussLongKokkos<DeviceType>::allocate()
{
  PairLJSwitch3CoulGaussLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_lj_sw3_cg**,Kokkos::LayoutRight,DeviceType>(
    "PairLJSwitch3CoulGaussLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSwitch3CoulGaussLongKokkos<DeviceType>::init_style()
{
  PairLJSwitch3CoulGaussLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJSwitch3CoulGaussLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJSwitch3CoulGaussLong::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i,j).lj2       = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3       = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4       = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).offset    = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cut_lj    = static_cast<KK_FLOAT>(cut_lj[i][j]);
  k_params.view_host()(i,j).cut_ljsq  = static_cast<KK_FLOAT>(cut_ljsqm);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_ljsqm;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsq;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cut_ljsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJSwitch3CoulGaussLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJSwitch3CoulGaussLongKokkos<LMPHostType>;
#endif
}
