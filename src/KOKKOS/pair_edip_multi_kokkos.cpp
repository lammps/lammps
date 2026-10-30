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

#include "pair_edip_multi_kokkos.h"

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

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairEDIPMultiKokkos<DeviceType>::PairEDIPMultiKokkos(LAMMPS *lmp) : PairEDIPMulti(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairEDIPMultiKokkos<DeviceType>::~PairEDIPMultiKokkos()
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
void PairEDIPMultiKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

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

  // build short neighbor list and the coordination of each atom

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum))
    d_neighbors_short = typename AT::t_int_2d_dl("EDIPMulti::neighbors_short",ignum*1.2,max_neighs);
  if ((int)d_numneigh_short.extent(0) < ignum) {
    d_numneigh_short = typename AT::t_int_1d("EDIPMulti::numneighs_short",ignum*1.2);
    d_zeta = typename AT::t_kkfloat_1d("EDIPMulti::zeta",ignum*1.2);
  }

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairEDIPMultiPreLoop>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairEDIPMultiPreLoop>(0,inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEDIPMultiCompute<HALFTHREAD,0> >(0,inum),*this);
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
   pre-loop: neighbors within the pair cutoff and environment coordination Z
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::operator()(TagPairEDIPMultiPreLoop, const int& ii) const {
  const int i = d_ilist[ii];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  KK_FLOAT zeta_i = 0.0;
  const int jnum = d_numneigh[i];
  int inside = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const int jtype = d_map[type[j]];
    const Param &param = d_params[d_elem3param(itype,jtype,jtype)];

    const KK_FLOAT delx = x(j,0) - xtmp;
    const KK_FLOAT dely = x(j,1) - ytmp;
    const KK_FLOAT delz = x(j,2) - ztmp;
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > static_cast<KK_FLOAT>(param.cutsq)) continue;

    d_neighbors_short(ii,inside) = j;
    inside++;

    const KK_FLOAT r_ij = Kokkos::sqrt(rsq);
    if (r_ij < static_cast<KK_FLOAT>(param.cutoffC))
      zeta_i += static_cast<KK_FLOAT>(1.0);
    else {
      KK_FLOAT fc, fdr;
      edip_fc(r_ij, param, fc, fdr);
      zeta_i += fc;
    }
  }
  d_numneigh_short(ii) = inside;
  d_zeta(ii) = zeta_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::operator()(TagPairEDIPMultiCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT zeta_i = d_zeta(ii);
  const int jnum = d_numneigh_short(ii);

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  KK_FLOAT dpairZ = 0.0;
  KK_FLOAT dtripleZ = 0.0;

  // two-body interactions
  // since pair i-j is different from pair j-i, double counting is
  // already considered in constructing the potential

  for (int jj = 0; jj < jnum; jj++) {
    KK_FLOAT dr_ij[3];
    KK_ACC_FLOAT f_ij[3];
    const int j = d_neighbors_short(ii,jj);
    const int jtype = d_map[type[j]];
    const Param &param_ij = d_params[d_elem3param(itype,jtype,jtype)];

    dr_ij[0] = x(j,0) - xtmp;
    dr_ij[1] = x(j,1) - ytmp;
    dr_ij[2] = x(j,2) - ztmp;
    const KK_FLOAT r_ij = Kokkos::sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2]);

    KK_FLOAT evdwl, fdr, fdZ;
    edip_pair(r_ij, zeta_i, param_ij, evdwl, fdr, fdZ);
    const KK_FLOAT fpair = -fdr / r_ij;
    dpairZ += fdZ;

    fxtmpi -= static_cast<KK_ACC_FLOAT>(fpair * dr_ij[0]);
    fytmpi -= static_cast<KK_ACC_FLOAT>(fpair * dr_ij[1]);
    fztmpi -= static_cast<KK_ACC_FLOAT>(fpair * dr_ij[2]);

    KK_ACC_FLOAT fxtmpj = static_cast<KK_ACC_FLOAT>(fpair * dr_ij[0]);
    KK_ACC_FLOAT fytmpj = static_cast<KK_ACC_FLOAT>(fpair * dr_ij[1]);
    KK_ACC_FLOAT fztmpj = static_cast<KK_ACC_FLOAT>(fpair * dr_ij[2]);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,-dr_ij[0],-dr_ij[1],-dr_ij[2]);
    }

    KK_FLOAT v1, v2;
    edip_fcut3(r_ij, param_ij, v1, v2);

    // three-body forces

    for (int kk = jj + 1; kk < jnum; kk++) {
      KK_FLOAT dr_ik[3];
      KK_ACC_FLOAT f_ik[3];
      int k = d_neighbors_short(ii,kk);
      const int ktype = d_map[type[k]];
      const Param &param_ik = d_params[d_elem3param(itype,ktype,ktype)];
      const Param &param_ijk = d_params[d_elem3param(itype,jtype,ktype)];

      dr_ik[0] = x(k,0) - xtmp;
      dr_ik[1] = x(k,1) - ytmp;
      dr_ik[2] = x(k,2) - ztmp;
      const KK_FLOAT r_ik = Kokkos::sqrt(dr_ik[0]*dr_ik[0] + dr_ik[1]*dr_ik[1] + dr_ik[2]*dr_ik[2]);

      const KK_FLOAT costheta = (dr_ij[0]*dr_ik[0] + dr_ij[1]*dr_ik[1] + dr_ij[2]*dr_ik[2]) / r_ij / r_ik;

      KK_FLOAT v3, v4, v5, v6, v7;
      edip_fcut3(r_ik, param_ik, v3, v4);
      edip_h(costheta, zeta_i, param_ijk, v5, v6, v7);

      // potential energy and forces

      const KK_FLOAT evdwl3 = v1 * v3 * v5;
      dtripleZ += v1 * v3 * v7;

      const KK_FLOAT dhl = v1 * v3 * v6;

      // derivatives of cos(theta) with respect to the positions of j and k

      const KK_FLOAT invrijrik = static_cast<KK_FLOAT>(1.0) / r_ij / r_ik;
      const KK_FLOAT ctrij = costheta / r_ij / r_ij;
      const KK_FLOAT ctrik = costheta / r_ik / r_ik;

      KK_FLOAT fpair_ij = -v2 * v3 * v5 / r_ij;
      KK_FLOAT fpair_ik = -v1 * v4 * v5 / r_ik;

      for (int n = 0; n < 3; n++) {
        const KK_FLOAT drj = invrijrik * dr_ik[n] - ctrij * dr_ij[n];
        const KK_FLOAT drk = invrijrik * dr_ij[n] - ctrik * dr_ik[n];
        f_ij[n] = static_cast<KK_ACC_FLOAT>(-dhl * drj + fpair_ij * dr_ij[n]);
        f_ik[n] = static_cast<KK_ACC_FLOAT>(-dhl * drk + fpair_ik * dr_ik[n]);
      }

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
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl3);
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl3,0.0,f_ij,f_ik,dr_ij,dr_ik);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  // forces due to environment coordination f(Z)

  const KK_FLOAT dZ = dpairZ + dtripleZ;
  for (int jj = 0; jj < jnum; jj++) {
    KK_FLOAT dr_ij[3];
    const int j = d_neighbors_short(ii,jj);
    const int jtype = d_map[type[j]];
    const Param &param = d_params[d_elem3param(itype,jtype,jtype)];

    dr_ij[0] = x(j,0) - xtmp;
    dr_ij[1] = x(j,1) - ytmp;
    dr_ij[2] = x(j,2) - ztmp;
    const KK_FLOAT r_ij = Kokkos::sqrt(dr_ij[0]*dr_ij[0] + dr_ij[1]*dr_ij[1] + dr_ij[2]*dr_ij[2]);
    if (r_ij < static_cast<KK_FLOAT>(param.cutoffC)) continue;

    KK_FLOAT fc, fdr;
    edip_fc(r_ij, param, fc, fdr);
    const KK_FLOAT dzetair = -fdr / r_ij * dZ;

    const KK_ACC_FLOAT fx = static_cast<KK_ACC_FLOAT>(dzetair * dr_ij[0]);
    const KK_ACC_FLOAT fy = static_cast<KK_ACC_FLOAT>(dzetair * dr_ij[1]);
    const KK_ACC_FLOAT fz = static_cast<KK_ACC_FLOAT>(dzetair * dr_ij[2]);

    a_f(j,0) += fx;
    a_f(j,1) += fy;
    a_f(j,2) += fz;

    fxtmpi -= fx;
    fytmpi -= fy;
    fztmpi -= fz;

    if (EVFLAG) {
      if (vflag_either)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,0.0,dzetair,-dr_ij[0],-dr_ij[1],-dr_ij[2]);
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
void PairEDIPMultiKokkos<DeviceType>::operator()(TagPairEDIPMultiCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairEDIPMultiCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   pair Vij, partial derivatives dVij(r,Z)/dr and dVij(r,Z)/dZ
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_pair(const KK_FLOAT &r, const KK_FLOAT &z,
                                                const Param &param, KK_FLOAT &eng,
                                                KK_FLOAT &fdr, KK_FLOAT &fZ) const
{
  const KK_FLOAT A = static_cast<KK_FLOAT>(param.A);
  const KK_FLOAT B = static_cast<KK_FLOAT>(param.B);
  const KK_FLOAT rho = static_cast<KK_FLOAT>(param.rho);
  const KK_FLOAT beta = static_cast<KK_FLOAT>(param.beta);
  KK_FLOAT v3, v4;

  const KK_FLOAT v1 = Kokkos::pow(B / r, rho);
  const KK_FLOAT v2 = Kokkos::exp(-beta * z * z);
  edip_fcut2(r, param, v3, v4);

  eng = A * (v1 - v2) * v3;
  fdr = A * (v1 - v2) * v4 + A * (-rho * v1 / r) * v3;
  fZ = A * (static_cast<KK_FLOAT>(2.0) * beta * z * v2) * v3;
}

/* ----------------------------------------------------------------------
   function fc(r) in calculating coordination Z and derivative fc'(r)
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_fc(const KK_FLOAT &r, const Param &param,
                                              KK_FLOAT &f_c, KK_FLOAT &fdr) const
{
  const KK_FLOAT a = static_cast<KK_FLOAT>(param.cutoffA);
  const KK_FLOAT c = static_cast<KK_FLOAT>(param.cutoffC);
  const KK_FLOAT alpha = static_cast<KK_FLOAT>(param.alpha);

  if (r < c + static_cast<KK_FLOAT>(1.0e-6)) {
    f_c = 1.0;
    fdr = 0.0;
    return;
  }

  if (r > a - static_cast<KK_FLOAT>(1.0e-6)) {
    f_c = 0.0;
    fdr = 0.0;
    return;
  }

  const KK_FLOAT xx = (a - c) / (r - c);
  const KK_FLOAT v1 = xx * xx * xx;
  const KK_FLOAT v2 = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) - v1);

  f_c = Kokkos::exp(alpha * v2);
  fdr = (static_cast<KK_FLOAT>(3.0) * xx * v1 / (a - c)) * (-alpha * v2 * v2) * f_c;
}

/* ----------------------------------------------------------------------
   cut-off function for Vij and its derivative fcut2'(r)
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_fcut2(const KK_FLOAT &r, const Param &param,
                                                 KK_FLOAT &f_c, KK_FLOAT &fdr) const
{
  const KK_FLOAT sigma = static_cast<KK_FLOAT>(param.sigma);
  const KK_FLOAT a = static_cast<KK_FLOAT>(param.cutoffA);

  if (r > a - static_cast<KK_FLOAT>(1.0e-6)) {
    f_c = 0.0;
    fdr = 0.0;
    return;
  }

  const KK_FLOAT v1 = static_cast<KK_FLOAT>(1.0) / (r - a);
  f_c = Kokkos::exp(sigma * v1);
  fdr = (-sigma * v1 * v1) * f_c;
}

/* ----------------------------------------------------------------------
   function tau(Z) and its derivative tau'(Z)
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_tau(const KK_FLOAT &z, const Param &param,
                                               KK_FLOAT &f_t, KK_FLOAT &fdZ) const
{
  const KK_FLOAT u1 = static_cast<KK_FLOAT>(param.u1);
  const KK_FLOAT u2 = static_cast<KK_FLOAT>(param.u2);
  const KK_FLOAT u3 = static_cast<KK_FLOAT>(param.u3);
  const KK_FLOAT u4 = static_cast<KK_FLOAT>(param.u4);

  const KK_FLOAT v1 = Kokkos::exp(-u4 * z);
  const KK_FLOAT v2 = Kokkos::exp(static_cast<KK_FLOAT>(-2.0) * u4 * z);

  f_t = u1 + u2 * u3 * v1 - u2 * v2;
  fdZ = -u2 * u3 * u4 * v1 + static_cast<KK_FLOAT>(2.0) * u2 * u4 * v2;
}

/* ----------------------------------------------------------------------
   function h(l,Z) and its partial derivatives dh(l,Z)/dl and dh(l,Z)/dZ
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_h(const KK_FLOAT &l, const KK_FLOAT &z,
                                             const Param &param, KK_FLOAT &f_h,
                                             KK_FLOAT &fdl, KK_FLOAT &fdZ) const
{
  const KK_FLOAT lambda = static_cast<KK_FLOAT>(param.lambda);
  const KK_FLOAT eta = static_cast<KK_FLOAT>(param.eta);
  const KK_FLOAT Q0 = static_cast<KK_FLOAT>(param.Q0);
  const KK_FLOAT mu = static_cast<KK_FLOAT>(param.mu);
  KK_FLOAT Tau, TaudZ;

  // function Q(Z) and its derivative Q'(Z)

  const KK_FLOAT Q = Q0 * Kokkos::exp(-mu * z);
  const KK_FLOAT QdZ = -mu * Q;

  edip_tau(z, param, Tau, TaudZ);

  const KK_FLOAT v1 = (l + Tau) * (l + Tau);
  const KK_FLOAT u2 = Q * v1;
  const KK_FLOAT v2 = Kokkos::exp(-u2);

  f_h = lambda * (static_cast<KK_FLOAT>(1.0) - v2 + eta * u2);

  // df/du2

  const KK_FLOAT v3 = lambda * (v2 + eta);

  // du2/dl

  const KK_FLOAT du2l = Q * static_cast<KK_FLOAT>(2.0) * (l + Tau);
  fdl = v3 * du2l;

  // du2/dZ

  const KK_FLOAT du2Z = QdZ * v1 + Q * static_cast<KK_FLOAT>(2.0) * (l + Tau) * TaudZ;
  fdZ = v3 * du2Z;
}

/* ----------------------------------------------------------------------
   cut-off function for Vijk and its derivative fcut3'(r)
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::edip_fcut3(const KK_FLOAT &r, const Param &param,
                                                 KK_FLOAT &f_c, KK_FLOAT &fdr) const
{
  const KK_FLOAT gamma = static_cast<KK_FLOAT>(param.gamma);
  const KK_FLOAT a = static_cast<KK_FLOAT>(param.cutoffA);

  if (r > a - static_cast<KK_FLOAT>(1.0e-6)) {
    f_c = 0.0;
    fdr = 0.0;
    return;
  }

  const KK_FLOAT v1 = static_cast<KK_FLOAT>(1.0) / (r - a);
  f_c = Kokkos::exp(gamma * v1);
  fdr = (-gamma * v1 * v1) * f_c;
}

/* ----------------------------------------------------------------------
   set coeffs and copy the parameters and element maps to the device
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEDIPMultiKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairEDIPMulti::coeff(narg,arg);

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.view_host();
  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];
  k_map.modify_host();
  k_map.template sync<DeviceType>();
  d_map = k_map.template view<DeviceType>();

  DAT::tdual_int_3d k_elem3param = DAT::tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  HAT::t_int_3d h_elem3param = k_elem3param.view_host();
  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem3param(i,j,k) = elem3param[i][j][k];
  k_elem3param.modify_host();
  k_elem3param.template sync<DeviceType>();
  d_elem3param = k_elem3param.template view<DeviceType>();

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.view_host();
  for (int m = 0; m < nparams; m++)
    h_params[m] = params[m];
  k_params.modify_host();
  k_params.template sync<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEDIPMultiKokkos<DeviceType>::init_style()
{
  PairEDIPMulti::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair edip/multi/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairEDIPMultiKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
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
void PairEDIPMultiKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
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
template class PairEDIPMultiKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairEDIPMultiKokkos<LMPHostType>;
#endif
}
