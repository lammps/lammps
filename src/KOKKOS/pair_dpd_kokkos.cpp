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
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#include "pair_dpd_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "memory_kokkos.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"
#include "atom_masks.h"
#include "kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double EPSILON = 1.0e-10;


template<class DeviceType>
PairDPDKokkos<DeviceType>::PairDPDKokkos(class LAMMPS *_lmp) :
  PairDPD(_lmp) ,
#ifdef DPD_USE_RAN_MARS
  rand_pool(0 /* unused */, _lmp)
#else
  rand_pool()
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDKokkos<DeviceType>::~PairDPDKokkos() {
  if (copymode) return;

#ifdef DPD_USE_RAN_MARS
  rand_pool.destroy();
#endif

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  memoryKK->destroy_kokkos(k_cutsq,cutsq);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDKokkos<DeviceType>::init_style()
{
  PairDPD::init_style();

#ifdef DPD_USE_RAN_MARS
  rand_pool.init(random,seed);
#else
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  rand_pool.init(seed + comm->me,unique_token.size());
#endif

  neighflag = lmp->kokkos->neighflag;

  if (force->newton_pair == 0 || neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style and newton on with pair dpd/kk");

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDKokkos<DeviceType>::compute(int eflagin, int vflagin)
{
  eflag = eflagin; vflag = vflagin;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,X_MASK | V_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();

  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_rf[0] = sqrt(force->special_lj[0]);
  special_rf[1] = sqrt(force->special_lj[1]);
  special_rf[2] = sqrt(force->special_lj[2]);
  special_rf[3] = sqrt(force->special_lj[3]);

  nlocal = atom->nlocal;
  dtinvsqrt = 1.0/sqrt(update->dt);

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

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

  // loop over neighbors of my atoms

  int inum = list->inum;
  EV_FLOAT ev;
  copymode = 1;
  if (neighflag == HALF) {
    if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDKokkos<HALF,1> >(0,inum),*this,ev);
    else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDKokkos<HALF,0> >(0,inum),*this);
  } else if (neighflag == HALFTHREAD) {
    if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDKokkos<HALFTHREAD,1> >(0,inum),*this,ev);
    else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDKokkos<HALFTHREAD,0> >(0,inum),*this);
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
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

  copymode = 0;

  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  else atomKK->modified(execution_space,F_MASK);

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_eatom = decltype(dup_eatom)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairDPDKokkos<DeviceType>::operator() (TagDPDKokkos<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagDPDKokkos<NEIGHFLAG,EVFLAG>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairDPDKokkos<DeviceType>::operator() (TagDPDKokkos<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT &ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int i,j,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,randnum,factor_dpd,factor_sqrt;
  double fx = 0,fy = 0,fz = 0;
  double evdwl = 0;
  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  vxtmp = v(i,0);
  vytmp = v(i,1);
  vztmp = v(i,2);
  itype = type(i);
  jnum = d_numneigh[i];
  rand_type rand_gen = rand_pool.get_state();
  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors(i,jj);
    factor_dpd = special_lj[sbmask(j)];
    factor_sqrt = special_rf[sbmask(j)];
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type(j);
    if (rsq < d_cutsq(itype,jtype)) {
      r = sqrt(rsq);
      if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
      rinv = 1.0/r;
      delvx = vxtmp - v(j,0);
      delvy = vytmp - v(j,1);
      delvz = vztmp - v(j,2);
      dot = delx*delvx + dely*delvy + delz*delvz;

      wd = 1.0 - r/params(itype,jtype).cut;

      randnum = rand_gen.normal();

      // conservative force
      fpair = params(itype,jtype).a0*wd;

      // drag force - parallel
      fpair -= params(itype,jtype).gamma*wd*wd*dot*rinv;
      fpair *= factor_dpd;

      // random force - parallel
      fpair += factor_sqrt*params(itype,jtype).sigma*wd*randnum*dtinvsqrt;
      fpair *= rinv;

      fx += fpair*delx;
      fy += fpair*dely;
      fz += fpair*delz;

      a_f(j,0) -= fpair*delx;
      a_f(j,1) -= fpair*dely;
      a_f(j,2) -= fpair*delz;

      if (EVFLAG && eflag_global) {
        // unshifted eng of conservative term:
        // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
        // eng shifted to 0.0 at cutoff
        evdwl = 0.5*params(itype,jtype).a0*params(itype,jtype).cut* wd*wd;
        evdwl *= factor_dpd;
        ev.evdwl += evdwl;
      }
      if (EVFLAG && (eflag_atom || vflag_either))
        this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }
  a_f(i,0) += fx;
  a_f(i,1) += fy;
  a_f(i,2) += fz;
  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairDPDKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

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
void PairDPDKokkos<DeviceType>::allocate()
{
  PairDPD::allocate();
  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType>("PairDPD::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairDPDKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDPDKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairDPD::init_one(i,j);

  k_params.h_view(i,j).cut = cut[i][j];
  k_params.h_view(i,j).a0 = a0[i][j];
  k_params.h_view(i,j).gamma = gamma[i][j];
  k_params.h_view(i,j).sigma = sigma[i][j];
  k_params.h_view(j,i) = k_params.h_view(i,j);

  k_params.template modify<LMPHostType>();

  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairDPDKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDPDKokkos<LMPHostType>;
#endif
}
