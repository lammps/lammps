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

#include "pair_buck6d_coul_gauss_long_kokkos.h"

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
PairBuck6dCoulGaussLongKokkos<DeviceType>::PairBuck6dCoulGaussLongKokkos(LAMMPS *lmp)
  : PairBuck6dCoulGaussLong(lmp)
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
PairBuck6dCoulGaussLongKokkos<DeviceType>::~PairBuck6dCoulGaussLongKokkos()
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
void PairBuck6dCoulGaussLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  EV_FLOAT ev;

  copymode = 1;

  ev = pair_compute<PairBuck6dCoulGaussLongKokkos<DeviceType>,void>
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
   Buck6d pair force with optional polynomial smoothing
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuck6dCoulGaussLongKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT r14inv = r6inv*r6inv*r2inv;

  const KK_FLOAT b1 = STACKPARAMS ? m_params[itype][jtype].buck6d1 : params(itype,jtype).buck6d1;
  const KK_FLOAT b2 = STACKPARAMS ? m_params[itype][jtype].buck6d2 : params(itype,jtype).buck6d2;
  const KK_FLOAT b3 = STACKPARAMS ? m_params[itype][jtype].buck6d3 : params(itype,jtype).buck6d3;
  const KK_FLOAT b4 = STACKPARAMS ? m_params[itype][jtype].buck6d4 : params(itype,jtype).buck6d4;

  const KK_FLOAT rexp = Kokkos::exp(-r*b2);
  const KK_FLOAT term1 = b3*r6inv;
  const KK_FLOAT term2 = b4*r14inv;
  const KK_FLOAT term3 = term2*term2;
  const KK_FLOAT term4 = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + term2);
  const KK_FLOAT term5 = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(2.0)*term2 + term3);
  KK_FLOAT forcebuck6d = b1*b2*r*rexp;
  forcebuck6d -= term1*(static_cast<KK_FLOAT>(6.0)*term4 - term5*static_cast<KK_FLOAT>(14.0)*term2);
  KK_FLOAT ebuck6d = b1*rexp - term1*term4;

  // optional polynomial smoothing near cutoff
  const KK_FLOAT rsmooth_sq_val = STACKPARAMS ? m_params[itype][jtype].rsmooth_sq : params(itype,jtype).rsmooth_sq;
  if (rsq > rsmooth_sq_val) {
    const KK_FLOAT rcu = r*rsq;
    const KK_FLOAT rqu = rsq*rsq;
    const KK_FLOAT c0v = STACKPARAMS ? m_params[itype][jtype].c0 : params(itype,jtype).c0;
    const KK_FLOAT c1v = STACKPARAMS ? m_params[itype][jtype].c1 : params(itype,jtype).c1;
    const KK_FLOAT c2v = STACKPARAMS ? m_params[itype][jtype].c2 : params(itype,jtype).c2;
    const KK_FLOAT c3v = STACKPARAMS ? m_params[itype][jtype].c3 : params(itype,jtype).c3;
    const KK_FLOAT c4v = STACKPARAMS ? m_params[itype][jtype].c4 : params(itype,jtype).c4;
    const KK_FLOAT c5v = STACKPARAMS ? m_params[itype][jtype].c5 : params(itype,jtype).c5;
    const KK_FLOAT sme = c5v*rqu*r + c4v*rqu + c3v*rcu + c2v*rsq + c1v*r + c0v;
    const KK_FLOAT smf = static_cast<KK_FLOAT>(5.0)*c5v*rqu + static_cast<KK_FLOAT>(4.0)*c4v*rcu
                       + static_cast<KK_FLOAT>(3.0)*c3v*rsq + static_cast<KK_FLOAT>(2.0)*c2v*r + c1v;
    forcebuck6d = forcebuck6d*sme - ebuck6d*smf*r;
  }

  return forcebuck6d*r2inv;
}

/* ----------------------------------------------------------------------
   Gaussian-damped Ewald Coulomb force with optional smoothing
   Uses polynomial erfc approximation in place of MathSpecial::my_erfcx
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuck6dCoulGaussLongKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;

  // Ewald factor: erf(g_ewald * r)
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT tg = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc_g = tg * (static_cast<KK_FLOAT>(A1)+tg*(static_cast<KK_FLOAT>(A2)+
                           tg * (static_cast<KK_FLOAT>(A3)+tg*(static_cast<KK_FLOAT>(A4)+
                           tg * static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT erf_ewald = static_cast<KK_FLOAT>(1.0) - erfc_g;

  // Gaussian-alpha factor: erf(alpha_ij * r)
  const KK_FLOAT alpha = STACKPARAMS ? m_params[itype][jtype].alpha_ij : params(itype,jtype).alpha_ij;
  const KK_FLOAT arg = alpha * r;
  const KK_FLOAT expa = Kokkos::exp(-arg*arg);
  const KK_FLOAT ta = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*arg);
  const KK_FLOAT erfc_a = ta * (static_cast<KK_FLOAT>(A1)+ta*(static_cast<KK_FLOAT>(A2)+
                           ta * (static_cast<KK_FLOAT>(A3)+ta*(static_cast<KK_FLOAT>(A4)+
                           ta * static_cast<KK_FLOAT>(A5))))) * expa;
  const KK_FLOAT erfa = static_cast<KK_FLOAT>(1.0) - erfc_a;

  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) * rinv;
  const KK_FLOAT falpha = erfa - static_cast<KK_FLOAT>(EWALD_F)*arg*expa;
  const KK_FLOAT ealpha = prefactor * (erfa - erf_ewald);

  KK_FLOAT forcecoul = prefactor * (falpha - erf_ewald + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor*falpha;

  // optional Coulomb smoothing near cutoff (global polynomial)
  if (rsq > rsmooth_sq_c_kk) {
    const KK_FLOAT rcu = r*rsq;
    const KK_FLOAT rqu = rsq*rsq;
    const KK_FLOAT sme = c5_c_kk*rqu*r + c4_c_kk*rqu + c3_c_kk*rcu
                       + c2_c_kk*rsq + c1_c_kk*r + c0_c_kk;
    const KK_FLOAT smf = static_cast<KK_FLOAT>(5.0)*c5_c_kk*rqu + static_cast<KK_FLOAT>(4.0)*c4_c_kk*rcu
                       + static_cast<KK_FLOAT>(3.0)*c3_c_kk*rsq + static_cast<KK_FLOAT>(2.0)*c2_c_kk*r + c1_c_kk;
    forcecoul = forcecoul*sme - ealpha*smf*r;
  }

  return forcecoul * r2inv;
}

/* ----------------------------------------------------------------------
   Buck6d pair energy with optional polynomial smoothing
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuck6dCoulGaussLongKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
               const int& itype, const int& jtype) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT r14inv = r6inv*r6inv*r2inv;

  const KK_FLOAT b1 = STACKPARAMS ? m_params[itype][jtype].buck6d1 : params(itype,jtype).buck6d1;
  const KK_FLOAT b2 = STACKPARAMS ? m_params[itype][jtype].buck6d2 : params(itype,jtype).buck6d2;
  const KK_FLOAT b3 = STACKPARAMS ? m_params[itype][jtype].buck6d3 : params(itype,jtype).buck6d3;
  const KK_FLOAT b4 = STACKPARAMS ? m_params[itype][jtype].buck6d4 : params(itype,jtype).buck6d4;
  const KK_FLOAT offset_val = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;

  const KK_FLOAT rexp = Kokkos::exp(-r*b2);
  const KK_FLOAT term1 = b3*r6inv;
  const KK_FLOAT term2 = b4*r14inv;
  const KK_FLOAT term4 = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + term2);
  KK_FLOAT ebuck6d = b1*rexp - term1*term4;

  const KK_FLOAT rsmooth_sq_val = STACKPARAMS ? m_params[itype][jtype].rsmooth_sq : params(itype,jtype).rsmooth_sq;
  if (rsq > rsmooth_sq_val) {
    const KK_FLOAT rcu = r*rsq;
    const KK_FLOAT rqu = rsq*rsq;
    const KK_FLOAT c0v = STACKPARAMS ? m_params[itype][jtype].c0 : params(itype,jtype).c0;
    const KK_FLOAT c1v = STACKPARAMS ? m_params[itype][jtype].c1 : params(itype,jtype).c1;
    const KK_FLOAT c2v = STACKPARAMS ? m_params[itype][jtype].c2 : params(itype,jtype).c2;
    const KK_FLOAT c3v = STACKPARAMS ? m_params[itype][jtype].c3 : params(itype,jtype).c3;
    const KK_FLOAT c4v = STACKPARAMS ? m_params[itype][jtype].c4 : params(itype,jtype).c4;
    const KK_FLOAT c5v = STACKPARAMS ? m_params[itype][jtype].c5 : params(itype,jtype).c5;
    const KK_FLOAT sme = c5v*rqu*r + c4v*rqu + c3v*rcu + c2v*rsq + c1v*r + c0v;
    ebuck6d *= sme;
  }

  return ebuck6d - offset_val;
}

/* ----------------------------------------------------------------------
   Gaussian-damped Ewald Coulomb energy with optional smoothing
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuck6dCoulGaussLongKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
               const int& itype, const int& jtype,
               const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;

  // Ewald factor: erf(g_ewald * r)
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT tg = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc_g = tg * (static_cast<KK_FLOAT>(A1)+tg*(static_cast<KK_FLOAT>(A2)+
                           tg * (static_cast<KK_FLOAT>(A3)+tg*(static_cast<KK_FLOAT>(A4)+
                           tg * static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT erf_ewald = static_cast<KK_FLOAT>(1.0) - erfc_g;

  // Gaussian-alpha factor: erf(alpha_ij * r)
  const KK_FLOAT alpha = STACKPARAMS ? m_params[itype][jtype].alpha_ij : params(itype,jtype).alpha_ij;
  const KK_FLOAT arg = alpha * r;
  const KK_FLOAT expa = Kokkos::exp(-arg*arg);
  const KK_FLOAT ta = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*arg);
  const KK_FLOAT erfc_a = ta * (static_cast<KK_FLOAT>(A1)+ta*(static_cast<KK_FLOAT>(A2)+
                           ta * (static_cast<KK_FLOAT>(A3)+ta*(static_cast<KK_FLOAT>(A4)+
                           ta * static_cast<KK_FLOAT>(A5))))) * expa;
  const KK_FLOAT erfa = static_cast<KK_FLOAT>(1.0) - erfc_a;

  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) * rinv;
  KK_FLOAT ealpha = prefactor * (erfa - erf_ewald);
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) ealpha -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor*erfa;

  // optional Coulomb smoothing near cutoff
  if (rsq > rsmooth_sq_c_kk) {
    const KK_FLOAT rcu = r*rsq;
    const KK_FLOAT rqu = rsq*rsq;
    const KK_FLOAT sme = c5_c_kk*rqu*r + c4_c_kk*rqu + c3_c_kk*rcu
                       + c2_c_kk*rsq + c1_c_kk*r + c0_c_kk;
    ealpha *= sme;
  }

  return ealpha;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuck6dCoulGaussLongKokkos<DeviceType>::allocate()
{
  PairBuck6dCoulGaussLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_buck6d**,Kokkos::LayoutRight,DeviceType>(
    "PairBuck6dCoulGaussLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuck6dCoulGaussLongKokkos<DeviceType>::init_style()
{
  PairBuck6dCoulGaussLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

  // copy global Coulomb smoothing parameters to KK_FLOAT class members
  c0_c_kk = static_cast<KK_FLOAT>(c0_c);
  c1_c_kk = static_cast<KK_FLOAT>(c1_c);
  c2_c_kk = static_cast<KK_FLOAT>(c2_c);
  c3_c_kk = static_cast<KK_FLOAT>(c3_c);
  c4_c_kk = static_cast<KK_FLOAT>(c4_c);
  c5_c_kk = static_cast<KK_FLOAT>(c5_c);
  rsmooth_sq_c_kk = static_cast<KK_FLOAT>(rsmooth_sq_c);

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
double PairBuck6dCoulGaussLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBuck6dCoulGaussLong::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i,j).buck6d1    = static_cast<KK_FLOAT>(buck6d1[i][j]);
  k_params.view_host()(i,j).buck6d2    = static_cast<KK_FLOAT>(buck6d2[i][j]);
  k_params.view_host()(i,j).buck6d3    = static_cast<KK_FLOAT>(buck6d3[i][j]);
  k_params.view_host()(i,j).buck6d4    = static_cast<KK_FLOAT>(buck6d4[i][j]);
  k_params.view_host()(i,j).c0         = static_cast<KK_FLOAT>(c0[i][j]);
  k_params.view_host()(i,j).c1         = static_cast<KK_FLOAT>(c1[i][j]);
  k_params.view_host()(i,j).c2         = static_cast<KK_FLOAT>(c2[i][j]);
  k_params.view_host()(i,j).c3         = static_cast<KK_FLOAT>(c3[i][j]);
  k_params.view_host()(i,j).c4         = static_cast<KK_FLOAT>(c4[i][j]);
  k_params.view_host()(i,j).c5         = static_cast<KK_FLOAT>(c5[i][j]);
  k_params.view_host()(i,j).rsmooth_sq = static_cast<KK_FLOAT>(rsmooth_sq[i][j]);
  k_params.view_host()(i,j).offset     = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).alpha_ij   = static_cast<KK_FLOAT>(alpha_ij[i][j]);
  k_params.view_host()(i,j).cut_ljsq   = static_cast<KK_FLOAT>(cut_ljsqm);
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
template class PairBuck6dCoulGaussLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBuck6dCoulGaussLongKokkos<LMPHostType>;
#endif
}
