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
------------------------------------------------------------------------- */

#include "pair_edip_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"

#include <cmath>
#include <limits>

using namespace LAMMPS_NS;
using namespace MathConst;

// same values as in pair_edip.cpp

static constexpr int GRIDDENSITY = 8000;
static constexpr double GRIDSTART = 0.1;
static constexpr int leadDimInteractionList = 64;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairEDIPKokkos<DeviceType>::PairEDIPKokkos(LAMMPS *lmp) : PairEDIP(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_overflow = DAT::tdual_int_scalar("pair:edip_overflow");
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairEDIPKokkos<DeviceType>::~PairEDIPKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEDIPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

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
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  cutsq_kk = static_cast<KK_FLOAT>(params[0].cutsq);

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  int chunk_size = 0;
  if (lmp->kokkos->threads_per_atom_set)
    chunk_size = lmp->kokkos->threads_per_atom;

  // build short neighbor list with the per-pair quantities and coordination of each atom

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = typename AT::t_int_2d_dl("EDIP::neighbors_short",ignum*1.2,max_neighs);
    d_preInvR = typename AT::t_kkfloat_2d_dl("EDIP::preInvR",ignum*1.2,max_neighs);
    d_preExp3B = typename AT::t_kkfloat_2d_dl("EDIP::preExp3B",ignum*1.2,max_neighs);
    d_preExp3BDerived = typename AT::t_kkfloat_2d_dl("EDIP::preExp3BDerived",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum) {
    d_numneigh_short = typename AT::t_int_1d("EDIP::numneighs_short",ignum*1.2);
    d_zeta = typename AT::t_kkfloat_1d("EDIP::zeta",ignum*1.2);
  }

  k_overflow.view_host()() = 0;
  k_overflow.modify_host();
  k_overflow.template sync<DeviceType>();
  d_overflow = k_overflow.template view<DeviceType>();

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairEDIPPreLoop>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairEDIPPreLoop>(0,inum), *this);

  // the interpolation grids of functions of the coordination have a fixed size

  k_overflow.template modify<DeviceType>();
  k_overflow.sync_host();
  if (k_overflow.view_host()())
    error->one(FLERR, "Pair style edip/kk supports at most {} neighbors within the cutoff per atom",
               leadDimInteractionList);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPCompute<HALFTHREAD,0> >(0,inum),*this);
    }
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += static_cast<double>(ev_all.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev_all.v[0]);
    virial[1] += static_cast<double>(ev_all.v[1]);
    virial[2] += static_cast<double>(ev_all.v[2]);
    virial[3] += static_cast<double>(ev_all.v[3]);
    virial[4] += static_cast<double>(ev_all.v[4]);
    virial[5] += static_cast<double>(ev_all.v[5]);
  }

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f            = {};
    dup_eatom        = {};
    dup_vatom        = {};
  }
}

/* ----------------------------------------------------------------------
   linear interpolation in a grid of the base class
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairEDIPKokkos<DeviceType>::interpolate(const typename AT::t_kkfloat_1d_randomread &grid,
                                                 const int &idx, const KK_FLOAT &frac) const
{
  const KK_FLOAT y1 = grid[idx];
  const KK_FLOAT y2 = grid[idx+1];
  return y1 + (y2 - y1) * frac;
}

/* ----------------------------------------------------------------------
   pre-loop: short neighbor list, per-pair quantities of the three-body
   term, and environment coordination Z
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPKokkos<DeviceType>::operator()(TagPairEDIPPreLoop, const int& ii) const {
  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  const KK_FLOAT cutoffA_kk = static_cast<KK_FLOAT>(cutoffA);
  const KK_FLOAT cutoffC_kk = static_cast<KK_FLOAT>(cutoffC);
  const KK_FLOAT gamm_kk = static_cast<KK_FLOAT>(gamm);

  KK_FLOAT zeta_i = 0.0;
  const int jnum = d_numneigh[i];
  int inside = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > cutsq_kk) continue;

    const KK_FLOAT r_ij = Kokkos::sqrt(rsq);
    const KK_FLOAT invR_ij = static_cast<KK_FLOAT>(1.0) / r_ij;
    const KK_FLOAT invRMinusCutoffA = static_cast<KK_FLOAT>(1.0) / (r_ij - cutoffA_kk);
    const KK_FLOAT gammInvRMinusCutoffA = gamm_kk * invRMinusCutoffA;

    const KK_FLOAT interpolTMP = (r_ij - static_cast<KK_FLOAT>(GRIDSTART)) * static_cast<KK_FLOAT>(GRIDDENSITY);
    const int interpolIDX = static_cast<int>(interpolTMP);
    const KK_FLOAT frac = interpolTMP - static_cast<KK_FLOAT>(interpolIDX);

    const KK_FLOAT exp3B_ij = interpolate(d_exp3B,interpolIDX,frac);

    d_neighbors_short(ii,inside) = j;
    d_preInvR(ii,inside) = invR_ij;
    d_preExp3B(ii,inside) = exp3B_ij;
    d_preExp3BDerived(ii,inside) = -exp3B_ij * gammInvRMinusCutoffA * invRMinusCutoffA;
    inside++;

    if (r_ij < cutoffC_kk)
      zeta_i += static_cast<KK_FLOAT>(1.0);
    else
      zeta_i += interpolate(d_cutoffFunction,interpolIDX,frac);
  }
  d_numneigh_short(ii) = inside;
  d_zeta(ii) = zeta_i;
  if (inside > leadDimInteractionList) d_overflow() = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPKokkos<DeviceType>::operator()(TagPairEDIPCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const KK_FLOAT cutoffA_kk = static_cast<KK_FLOAT>(cutoffA);
  const KK_FLOAT cutoffC_kk = static_cast<KK_FLOAT>(cutoffC);
  const KK_FLOAT sigma_kk = static_cast<KK_FLOAT>(sigma);
  const KK_FLOAT rho_kk = static_cast<KK_FLOAT>(rho);
  const KK_FLOAT lambda_kk = static_cast<KK_FLOAT>(lambda);
  const KK_FLOAT eta_kk = static_cast<KK_FLOAT>(eta);
  const KK_FLOAT mu_kk = static_cast<KK_FLOAT>(mu);
  const KK_FLOAT beta_kk = static_cast<KK_FLOAT>(beta);

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh_short[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  // quantities depending on zeta_i

  const KK_FLOAT zeta_i = d_zeta[ii];
  KK_FLOAT interpolTMP = zeta_i * static_cast<KK_FLOAT>(GRIDDENSITY);
  int interpolIDX = static_cast<int>(interpolTMP);
  KK_FLOAT frac = interpolTMP - static_cast<KK_FLOAT>(interpolIDX);

  const KK_FLOAT expMinusBetaZeta_iZeta_i = interpolate(d_expMinusBetaZeta_iZeta_iGrid,interpolIDX,frac);
  const KK_FLOAT qFunction = interpolate(d_qFunctionGrid,interpolIDX,frac);
  const KK_FLOAT tauFunction = interpolate(d_tauFunctionGrid,interpolIDX,frac);
  const KK_FLOAT tauFunctionDerived = interpolate(d_tauFunctionDerivedGrid,interpolIDX,frac);

  const KK_FLOAT forceModCoord_factor = static_cast<KK_FLOAT>(2.0) * beta_kk * zeta_i * expMinusBetaZeta_iZeta_i;

  KK_FLOAT forceModCoord = 0.0;

  // two-body and three-body interactions

  for (int jj = 0; jj < jnum; jj++) {
    KK_FLOAT dr_ij[3];
    KK_ACC_FLOAT f_ij[3];
    const int j = d_neighbors_short(ii,jj);

    dr_ij[0] = x(j,0) - xtmp;
    dr_ij[1] = x(j,1) - ytmp;
    dr_ij[2] = x(j,2) - ztmp;
    const KK_FLOAT r_ij = Kokkos::sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2]);

    // two-body quantities, interpolated again instead of stored

    interpolTMP = (r_ij - static_cast<KK_FLOAT>(GRIDSTART)) * static_cast<KK_FLOAT>(GRIDDENSITY);
    interpolIDX = static_cast<int>(interpolTMP);
    frac = interpolTMP - static_cast<KK_FLOAT>(interpolIDX);

    const KK_FLOAT invRMinusCutoffA = static_cast<KK_FLOAT>(1.0) / (r_ij - cutoffA_kk);
    const KK_FLOAT sigmaInvRMinusCutoffA = sigma_kk * invRMinusCutoffA;
    const KK_FLOAT exp2B_ij = interpolate(d_exp2B,interpolIDX,frac);
    const KK_FLOAT exp2BDerived_ij = -exp2B_ij * sigmaInvRMinusCutoffA * invRMinusCutoffA;
    const KK_FLOAT pow2B_ij = interpolate(d_pow2B,interpolIDX,frac);

    const KK_FLOAT invR_ij = d_preInvR(ii,jj);

    const KK_FLOAT potential2B_factor = pow2B_ij - expMinusBetaZeta_iZeta_i;
    const KK_FLOAT pow2BDerived_ij = -rho_kk * invR_ij * pow2B_ij;

    forceModCoord += (forceModCoord_factor * exp2B_ij);

    const KK_FLOAT forceMod2B = exp2BDerived_ij * potential2B_factor + exp2B_ij * pow2BDerived_ij;

    const KK_FLOAT directorCos_ij_x = invR_ij * dr_ij[0];
    const KK_FLOAT directorCos_ij_y = invR_ij * dr_ij[1];
    const KK_FLOAT directorCos_ij_z = invR_ij * dr_ij[2];

    const KK_FLOAT exp3B_ij = d_preExp3B(ii,jj);
    const KK_FLOAT exp3BDerived_ij = d_preExp3BDerived(ii,jj);

    f_ij[0] = static_cast<KK_ACC_FLOAT>(forceMod2B * directorCos_ij_x);
    f_ij[1] = static_cast<KK_ACC_FLOAT>(forceMod2B * directorCos_ij_y);
    f_ij[2] = static_cast<KK_ACC_FLOAT>(forceMod2B * directorCos_ij_z);

    fxtmpi += f_ij[0];
    fytmpi += f_ij[1];
    fztmpi += f_ij[2];

    KK_ACC_FLOAT fxtmpj = -f_ij[0];
    KK_ACC_FLOAT fytmpj = -f_ij[1];
    KK_ACC_FLOAT fztmpj = -f_ij[2];

    if (EVFLAG) {
      const KK_FLOAT evdwl = exp2B_ij * potential2B_factor;
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,-forceMod2B*invR_ij,dr_ij[0],dr_ij[1],dr_ij[2]);
    }

    // three-body forces

    for (int kk = jj + 1; kk < jnum; kk++) {
      KK_FLOAT dr_ik[3];
      KK_ACC_FLOAT f_ik[3];
      int k = d_neighbors_short(ii,kk);

      dr_ik[0] = x(k,0) - xtmp;
      dr_ik[1] = x(k,1) - ytmp;
      dr_ik[2] = x(k,2) - ztmp;

      const KK_FLOAT invR_ik = d_preInvR(ii,kk);

      const KK_FLOAT directorCos_ik_x = invR_ik * dr_ik[0];
      const KK_FLOAT directorCos_ik_y = invR_ik * dr_ik[1];
      const KK_FLOAT directorCos_ik_z = invR_ik * dr_ik[2];

      const KK_FLOAT cosTeta = directorCos_ij_x * directorCos_ik_x + directorCos_ij_y * directorCos_ik_y +
          directorCos_ij_z * directorCos_ik_z;

      const KK_FLOAT cosTetaDiff = cosTeta + tauFunction;
      const KK_FLOAT cosTetaDiffCosTetaDiff = cosTetaDiff * cosTetaDiff;
      const KK_FLOAT qFunctionCosTetaDiffCosTetaDiff = cosTetaDiffCosTetaDiff * qFunction;
      const KK_FLOAT expMinusQFunctionCosTetaDiffCosTetaDiff = Kokkos::exp(-qFunctionCosTetaDiffCosTetaDiff);

      const KK_FLOAT potentia3B_factor = lambda_kk *
          ((static_cast<KK_FLOAT>(1.0) - expMinusQFunctionCosTetaDiffCosTetaDiff) +
           eta_kk * qFunctionCosTetaDiffCosTetaDiff);

      const KK_FLOAT exp3B_ik = d_preExp3B(ii,kk);
      const KK_FLOAT exp3BDerived_ik = d_preExp3BDerived(ii,kk);

      const KK_FLOAT forceMod3B_factor1_ij = -exp3BDerived_ij * exp3B_ik * potentia3B_factor;
      const KK_FLOAT forceMod3B_factor2 = static_cast<KK_FLOAT>(2.0) * lambda_kk * exp3B_ij * exp3B_ik *
          qFunction * cosTetaDiff * (eta_kk + expMinusQFunctionCosTetaDiffCosTetaDiff);
      const KK_FLOAT forceMod3B_factor2_ij = forceMod3B_factor2 * invR_ij;

      f_ij[0] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ij * directorCos_ij_x +
          forceMod3B_factor2_ij * (cosTeta * directorCos_ij_x - directorCos_ik_x));
      f_ij[1] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ij * directorCos_ij_y +
          forceMod3B_factor2_ij * (cosTeta * directorCos_ij_y - directorCos_ik_y));
      f_ij[2] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ij * directorCos_ij_z +
          forceMod3B_factor2_ij * (cosTeta * directorCos_ij_z - directorCos_ik_z));

      const KK_FLOAT forceMod3B_factor1_ik = -exp3BDerived_ik * exp3B_ij * potentia3B_factor;
      const KK_FLOAT forceMod3B_factor2_ik = forceMod3B_factor2 * invR_ik;

      f_ik[0] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ik * directorCos_ik_x +
          forceMod3B_factor2_ik * (cosTeta * directorCos_ik_x - directorCos_ij_x));
      f_ik[1] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ik * directorCos_ik_y +
          forceMod3B_factor2_ik * (cosTeta * directorCos_ik_y - directorCos_ij_y));
      f_ik[2] = static_cast<KK_ACC_FLOAT>(forceMod3B_factor1_ik * directorCos_ik_z +
          forceMod3B_factor2_ik * (cosTeta * directorCos_ik_z - directorCos_ij_z));

      forceModCoord += (forceMod3B_factor2 * (tauFunctionDerived - static_cast<KK_FLOAT>(0.5) * mu_kk * cosTetaDiff));

      fxtmpj += f_ij[0];
      fytmpj += f_ij[1];
      fztmpj += f_ij[2];

      a_f(k,0) += f_ik[0];
      a_f(k,1) += f_ik[1];
      a_f(k,2) += f_ik[2];

      fxtmpi -= f_ij[0] + f_ik[0];
      fytmpi -= f_ij[1] + f_ik[1];
      fztmpi -= f_ij[2] + f_ik[2];

      if (EVFLAG) {
        const KK_FLOAT evdwl = exp3B_ij * exp3B_ik * potentia3B_factor;
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,f_ij,f_ik,dr_ij,dr_ik);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  // forces due to environment coordination f(Z)

  for (int jj = 0; jj < jnum; jj++) {
    KK_FLOAT dr_ij[3];
    KK_ACC_FLOAT f_ij[3];
    const int j = d_neighbors_short(ii,jj);

    dr_ij[0] = xtmp - x(j,0);
    dr_ij[1] = ytmp - x(j,1);
    dr_ij[2] = ztmp - x(j,2);
    const KK_FLOAT r_ij = Kokkos::sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2]);
    if (r_ij < cutoffC_kk) continue;

    interpolTMP = (r_ij - static_cast<KK_FLOAT>(GRIDSTART)) * static_cast<KK_FLOAT>(GRIDDENSITY);
    interpolIDX = static_cast<int>(interpolTMP);
    frac = interpolTMP - static_cast<KK_FLOAT>(interpolIDX);

    const KK_FLOAT zeta_iDerived = interpolate(d_cutoffFunctionDerived,interpolIDX,frac);
    const KK_FLOAT zeta_iDerivedInvR_ij = zeta_iDerived * d_preInvR(ii,jj);
    const KK_FLOAT forceModCoord_ij = forceModCoord * zeta_iDerivedInvR_ij;

    f_ij[0] = static_cast<KK_ACC_FLOAT>(forceModCoord_ij * dr_ij[0]);
    f_ij[1] = static_cast<KK_ACC_FLOAT>(forceModCoord_ij * dr_ij[1]);
    f_ij[2] = static_cast<KK_ACC_FLOAT>(forceModCoord_ij * dr_ij[2]);

    fxtmpi -= f_ij[0];
    fytmpi -= f_ij[1];
    fztmpi -= f_ij[2];

    a_f(j,0) += f_ij[0];
    a_f(j,1) += f_ij[1];
    a_f(j,2) += f_ij[2];

    if (EVFLAG) {
      if (vflag_either)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,0.0,-forceModCoord_ij,dr_ij[0],dr_ij[1],dr_ij[2]);
    }
  }

  a_f(i,0) += fxtmpi;
  a_f(i,1) += fytmpi;
  a_f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPKokkos<DeviceType>::operator()(TagPairEDIPCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairEDIPCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs and copy the interpolation grids to the device
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEDIPKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairEDIP::coeff(narg,arg);

  // grid sizes as in PairEDIP::allocateGrids()

  double const leftLimitToZero = -std::numeric_limits<double>::min() * 1000.0;
  const int numGridPointsZ = (int) (leadDimInteractionList * GRIDDENSITY) + 2;
  const int numGridPointsCutoffFunction = (int) ((cutoffC - GRIDSTART) * GRIDDENSITY) +
    (int) ((cutoffA - cutoffC) * GRIDDENSITY) + 2;
  const int numGridPointsRTotal = (int) ((cutoffA + leftLimitToZero - GRIDSTART) * GRIDDENSITY) + 2;

  struct {
    DAT::tdual_kkfloat_1d &kview;
    typename AT::t_kkfloat_1d_randomread &dview;
    const double *grid;
    int size;
    const char *name;
  } grids[] = {
    {k_cutoffFunction, d_cutoffFunction, cutoffFunction, numGridPointsCutoffFunction, "edip:cutoffFunction"},
    {k_cutoffFunctionDerived, d_cutoffFunctionDerived, cutoffFunctionDerived, numGridPointsCutoffFunction,
     "edip:cutoffFunctionDerived"},
    {k_pow2B, d_pow2B, pow2B, numGridPointsRTotal, "edip:pow2B"},
    {k_exp2B, d_exp2B, exp2B, numGridPointsRTotal, "edip:exp2B"},
    {k_exp3B, d_exp3B, exp3B, numGridPointsRTotal, "edip:exp3B"},
    {k_qFunctionGrid, d_qFunctionGrid, qFunctionGrid, numGridPointsZ, "edip:qFunctionGrid"},
    {k_expMinusBetaZeta_iZeta_iGrid, d_expMinusBetaZeta_iZeta_iGrid, expMinusBetaZeta_iZeta_iGrid,
     numGridPointsZ, "edip:expMinusBetaZeta_iZeta_iGrid"},
    {k_tauFunctionGrid, d_tauFunctionGrid, tauFunctionGrid, numGridPointsZ, "edip:tauFunctionGrid"},
    {k_tauFunctionDerivedGrid, d_tauFunctionDerivedGrid, tauFunctionDerivedGrid, numGridPointsZ,
     "edip:tauFunctionDerivedGrid"}
  };

  for (auto &g : grids) {
    g.kview = DAT::tdual_kkfloat_1d(g.name,g.size);
    auto h_view = g.kview.view_host();
    for (int n = 0; n < g.size; n++) h_view(n) = static_cast<KK_FLOAT>(g.grid[n]);
    g.kview.modify_host();
    g.kview.template sync<DeviceType>();
    g.dview = g.kview.template view<DeviceType>();
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEDIPKokkos<DeviceType>::init_style()
{
  PairEDIP::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair edip/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    const KK_FLOAT epairhalf = static_cast<KK_FLOAT>(0.5) * epair;
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairhalf);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairhalf);
  }

  if (vflag_either) {
    KK_FLOAT v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (vflag_global) {
      for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(v[n]);
    }

    if (vflag_atom) {
      for (int n = 0; n < 6; n++) {
        const KK_ACC_FLOAT vhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v[n]);
        a_vatom(i,n) += vhalf;
        a_vatom(j,n) += vhalf;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
          const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                     KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const
{
  KK_FLOAT epairthird,v[6];

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    epairthird = static_cast<KK_FLOAT>(THIRD) * (evdwl + ecoul);
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[k] += static_cast<KK_ACC_FLOAT>(epairthird);
  }

  if (vflag_either) {
    v[0] = drji[0]*static_cast<KK_FLOAT>(fj[0]) + drki[0]*static_cast<KK_FLOAT>(fk[0]);
    v[1] = drji[1]*static_cast<KK_FLOAT>(fj[1]) + drki[1]*static_cast<KK_FLOAT>(fk[1]);
    v[2] = drji[2]*static_cast<KK_FLOAT>(fj[2]) + drki[2]*static_cast<KK_FLOAT>(fk[2]);
    v[3] = drji[0]*static_cast<KK_FLOAT>(fj[1]) + drki[0]*static_cast<KK_FLOAT>(fk[1]);
    v[4] = drji[0]*static_cast<KK_FLOAT>(fj[2]) + drki[0]*static_cast<KK_FLOAT>(fk[2]);
    v[5] = drji[1]*static_cast<KK_FLOAT>(fj[2]) + drki[1]*static_cast<KK_FLOAT>(fk[2]);

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v[0]);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v[1]);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v[2]);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v[3]);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v[4]);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v[5]);
    }

    if (vflag_atom) {
      for (int n = 0; n < 6; n++) {
        const KK_ACC_FLOAT vthird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[n]);
        a_vatom(i,n) += vthird;
        a_vatom(j,n) += vthird;
        a_vatom(k,n) += vthird;
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairEDIPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairEDIPKokkos<LMPHostType>;
#endif
}
