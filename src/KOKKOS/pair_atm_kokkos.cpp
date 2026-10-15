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
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_atm_kokkos.h"

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
PairATMKokkos<DeviceType>::PairATMKokkos(LAMMPS *lmp) : PairATM(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  nu_synced = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairATMKokkos<DeviceType>::~PairATMKokkos()
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
void PairATMKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  // copy the type triplet coefficients to the device after init_one() changed them

  if (!nu_synced) {
    const int np1 = atom->ntypes + 1;
    k_nu = DAT::tdual_kkfloat_3d("pair:nu",np1,np1,np1);
    auto h_nu = k_nu.view_host();
    for (int i = 1; i < np1; i++)
      for (int j = 1; j < np1; j++)
        for (int k = 1; k < np1; k++)
          h_nu(i,j,k) = static_cast<KK_FLOAT>(nu[i][j][k]);
    k_nu.modify_host();
    k_nu.template sync<DeviceType>();
    d_nu = k_nu.template view<DeviceType>();
    nu_synced = 1;
  }

  cutoff_squared = static_cast<KK_FLOAT>(cut_global*cut_global);
  const double triple = cut_triple*cut_triple*cut_triple;
  cutoff_triple_sixth = static_cast<KK_FLOAT>(triple*triple);

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

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = typename AT::t_int_2d_dl("ATM::neighbors_short",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum)
    d_numneigh_short = typename AT::t_int_1d("ATM::numneighs_short",ignum*1.2);

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairATMComputeShortNeigh>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairATMComputeShortNeigh>(0,inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairATMCompute<HALFTHREAD,0> >(0,inum),*this);
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
   short neighbor list: neighbors within the cutoff that are in the half
   space of larger x (ties broken by y and z), so that each I-J-K triplet
   is computed once, by the owner of the atom with the smallest x coord
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairATMKokkos<DeviceType>::operator()(TagPairATMComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const KK_FLOAT xi = x(i,0);
    const KK_FLOAT yi = x(i,1);
    const KK_FLOAT zi = x(i,2);

    const int jnum = d_numneigh[i];
    int inside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const KK_FLOAT delx = x(j,0) - xi;
      if (delx < static_cast<KK_FLOAT>(0.0)) continue;
      const KK_FLOAT dely = x(j,1) - yi;
      if (delx == static_cast<KK_FLOAT>(0.0) && dely < static_cast<KK_FLOAT>(0.0)) continue;
      const KK_FLOAT delz = x(j,2) - zi;
      if (delx == static_cast<KK_FLOAT>(0.0) && dely == static_cast<KK_FLOAT>(0.0) &&
          delz < static_cast<KK_FLOAT>(0.0)) continue;
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutoff_squared) {
        d_neighbors_short(ii,inside) = j;
        inside++;
      }
    }
    d_numneigh_short(ii) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairATMKokkos<DeviceType>::operator()(TagPairATMCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  KK_FLOAT rij[3],rik[3],rjk[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;

  const int i = d_ilist[ii];
  const int itype = type[i];
  const KK_FLOAT xi = x(i,0);
  const KK_FLOAT yi = x(i,1);
  const KK_FLOAT zi = x(i,2);

  const int jnum = d_numneigh_short[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  const int jnumm1 = jnum - 1;

  for (int jj = 0; jj < jnumm1; jj++) {
    const int j = d_neighbors_short(ii,jj);
    const int jtype = type[j];
    rij[0] = x(j,0) - xi;
    rij[1] = x(j,1) - yi;
    rij[2] = x(j,2) - zi;
    const KK_FLOAT rij2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

    KK_ACC_FLOAT fxtmpj = 0.0;
    KK_ACC_FLOAT fytmpj = 0.0;
    KK_ACC_FLOAT fztmpj = 0.0;

    for (int kk = jj+1; kk < jnum; kk++) {
      int k = d_neighbors_short(ii,kk);

      rik[0] = x(k,0) - xi;
      rik[1] = x(k,1) - yi;
      rik[2] = x(k,2) - zi;
      const KK_FLOAT rik2 = rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2];

      rjk[0] = x(k,0) - x(j,0);
      rjk[1] = x(k,1) - x(j,1);
      rjk[2] = x(k,2) - x(j,2);
      const KK_FLOAT rjk2 = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2];
      if (rjk2 > cutoff_squared) continue;

      const KK_FLOAT r6 = rij2*rjk2*rik2;
      if (r6 > cutoff_triple_sixth) continue;

      const KK_FLOAT nu_local = d_nu(itype,jtype,type[k]);
      if (nu_local == static_cast<KK_FLOAT>(0.0)) continue;

      interaction_ddd_kk(nu_local,r6,rij2,rik2,rjk2,rij,rik,rjk,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];
      fxtmpj += fj[0];
      fytmpj += fj[1];
      fztmpj += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,fj,fk,rij,rik);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  a_f(i,0) += fxtmpi;
  a_f(i,1) += fytmpi;
  a_f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairATMKokkos<DeviceType>::operator()(TagPairATMCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairATMCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairATMKokkos<DeviceType>::init_style()
{
  PairATM::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair atm/kk");
}

/* ----------------------------------------------------------------------
   init for one i,j type pair; the base class sets all permutations of
   nu, so the device copy must be refreshed
------------------------------------------------------------------------- */

template<class DeviceType>
double PairATMKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairATM::init_one(i,j);
  nu_synced = 0;
  return cutone;
}

/* ----------------------------------------------------------------------
   Axilrod-Teller-Muto (dipole-dipole-dipole) potential
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairATMKokkos<DeviceType>::interaction_ddd_kk(const KK_FLOAT &nu, const KK_FLOAT &r6,
                              const KK_FLOAT &rij2, const KK_FLOAT &rik2, const KK_FLOAT &rjk2,
                              KK_FLOAT *rij, KK_FLOAT *rik, KK_FLOAT *rjk,
                              KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, const int &eflag, KK_FLOAT &eng) const
{
  const KK_FLOAT r5inv = nu / (r6*r6*Kokkos::sqrt(r6));
  const KK_FLOAT rri = rik[0]*rij[0] + rik[1]*rij[1] + rik[2]*rij[2];
  const KK_FLOAT rrj = rij[0]*rjk[0] + rij[1]*rjk[1] + rij[2]*rjk[2];
  const KK_FLOAT rrk = rjk[0]*rik[0] + rjk[1]*rik[1] + rjk[2]*rik[2];
  const KK_FLOAT rrr = static_cast<KK_FLOAT>(5.0)*rri*rrj*rrk;
  const KK_FLOAT fac = static_cast<KK_FLOAT>(3.0)*r5inv;
  for (int n = 0; n < 3; n++) {
    const KK_FLOAT fjn = rrj*(rrk - rri)*rik[n] -
      (rrk*rri - rjk2*rik2 + rrr/rij2) * rij[n] +
      (rrk*rri - rik2*rij2 + rrr/rjk2) * rjk[n];
    fj[n] = static_cast<KK_ACC_FLOAT>(fjn*fac);
    const KK_FLOAT fkn = rrk*(rri + rrj)*rij[n] +
      (rri*rrj + rik2*rij2 - rrr/rjk2) * rjk[n] +
      (rri*rrj + rij2*rjk2 - rrr/rik2) * rik[n];
    fk[n] = static_cast<KK_ACC_FLOAT>(fkn*fac);
  }
  if (eflag) eng = (r6 - static_cast<KK_FLOAT>(0.6)*rrr)*r5inv;
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
void PairATMKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
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
template class PairATMKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairATMKokkos<LMPHostType>;
#endif
}
