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

#include "pair_dpd_coul_slater_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDCoulSlaterLongKokkos<DeviceType>::PairDPDCoulSlaterLongKokkos(class LAMMPS *_lmp) :
  PairDPDCoulSlaterLong(_lmp) ,
#ifdef DPD_USE_RAN_MARS
  rand_pool(0 /* unused */, _lmp)
#else
  rand_pool()
#endif
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDCoulSlaterLongKokkos<DeviceType>::~PairDPDCoulSlaterLongKokkos() {
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
void PairDPDCoulSlaterLongKokkos<DeviceType>::init_style()
{
  PairDPDCoulSlaterLong::init_style();

#ifdef DPD_USE_RAN_MARS
  rand_pool.init(random,seed);
#else
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  rand_pool = decltype(rand_pool)(seed + comm->me,unique_token.size());
#endif

  neighflag = lmp->kokkos->neighflag;

  if (force->newton_pair == 0 || neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style and newton on with pair dpd/coul/slater/long/kk");

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDCoulSlaterLongKokkos<DeviceType>::compute(int eflagin, int vflagin)
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

  atomKK->sync(execution_space,X_MASK | V_MASK | F_MASK | Q_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();

  for (int m = 0; m < 4; ++m) {
    special_lj[m] = static_cast<KK_FLOAT>(force->special_lj[m]);
    special_rf[m] = static_cast<KK_FLOAT>(sqrt(force->special_lj[m]));
    special_coul[m] = static_cast<KK_FLOAT>(force->special_coul[m]);
  }

  nlocal = atom->nlocal;
  dtinvsqrt = static_cast<KK_FLOAT>(1.0/sqrt(update->dt));
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  lamda_kk = static_cast<KK_FLOAT>(lamda);

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
    if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDCoulSlaterLongKokkos<HALF,1> >(0,inum),*this,ev);
    else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDCoulSlaterLongKokkos<HALF,0> >(0,inum),*this);
  } else if (neighflag == HALFTHREAD) {
    if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDCoulSlaterLongKokkos<HALFTHREAD,1> >(0,inum),*this,ev);
    else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDCoulSlaterLongKokkos<HALFTHREAD,0> >(0,inum),*this);
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) {
    eng_vdwl += static_cast<double>(ev.evdwl);
    eng_coul += static_cast<double>(ev.ecoul);
  }
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

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

  copymode = 0;

  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  else atomKK->modified(execution_space,F_MASK);

  // free duplicated memory
  if (need_dup) {
    dup_f     = {};
    dup_eatom = {};
    dup_vatom = {};
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairDPDCoulSlaterLongKokkos<DeviceType>::operator() (TagDPDCoulSlaterLongKokkos<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagDPDCoulSlaterLongKokkos<NEIGHFLAG,EVFLAG>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairDPDCoulSlaterLongKokkos<DeviceType>::operator() (TagDPDCoulSlaterLongKokkos<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT &ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT vxtmp = v(i,0);
  const KK_FLOAT vytmp = v(i,1);
  const KK_FLOAT vztmp = v(i,2);
  const KK_FLOAT qtmp = q(i);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  KK_FLOAT fx = 0, fy = 0, fz = 0;

  rand_type rand_gen = rand_pool.get_state();
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_dpd = special_lj[sbmask(j)];
    const KK_FLOAT factor_sqrt = special_rf[sbmask(j)];
    const KK_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const int jtype = type(j);

    if (rsq < d_cutsq(itype,jtype)) {
      const KK_FLOAT r = Kokkos::sqrt(rsq);
      KK_FLOAT evdwl = 0;
      KK_FLOAT ecoul = 0;

      // DPD force if distance below DPD cutoff; r can be 0.0 in DPD systems

      KK_FLOAT forcedpd = 0;
      if ((rsq < params(itype,jtype).cut_dpdsq) && (r > static_cast<KK_FLOAT>(EPSILON))) {
        const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0)/r;
        const KK_FLOAT delvx = vxtmp - v(j,0);
        const KK_FLOAT delvy = vytmp - v(j,1);
        const KK_FLOAT delvz = vztmp - v(j,2);
        const KK_FLOAT dot = delx*delvx + dely*delvy + delz*delvz;
        const KK_FLOAT wd = static_cast<KK_FLOAT>(1.0) - r/params(itype,jtype).cut_dpd;
        const KK_FLOAT randnum = static_cast<KK_FLOAT>(rand_gen.normal());

        // conservative force = a0 * wd
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;
        // random force must be scaled by sqrt(factor_dpd)

        forcedpd = params(itype,jtype).a0*wd;
        forcedpd -= params(itype,jtype).gamma*wd*wd*dot*rinv;
        forcedpd *= factor_dpd;
        forcedpd += factor_sqrt*params(itype,jtype).sigma*wd*randnum*dtinvsqrt;
        forcedpd *= rinv;

        // energy shifted to 0.0 at the cutoff
        if (EVFLAG && eflag)
          evdwl = static_cast<KK_FLOAT>(0.5)*params(itype,jtype).a0*params(itype,jtype).cut_dpd*wd*wd*factor_dpd;
      }

      // Slater-smeared long-range Coulomb force if distance below Slater cutoff

      KK_FLOAT forcecoul = 0;
      if (rsq < params(itype,jtype).cut_slatersq) {
        const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
        const KK_FLOAT grij = g_ewald_kk * r;
        const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
        const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
        const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                              t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                              t*static_cast<KK_FLOAT>(A5))))) * expm2;
        const KK_FLOAT rlamda = r/lamda_kk;
        const KK_FLOAT exp2rl = Kokkos::exp(static_cast<KK_FLOAT>(-2.0)*rlamda);
        const KK_FLOAT slater_term = exp2rl*(static_cast<KK_FLOAT>(1.0) + (static_cast<KK_FLOAT>(2.0)*rlamda*(static_cast<KK_FLOAT>(1.0)+rlamda)));
        const KK_FLOAT prefactor = qqrd2e * qtmp*q(j)/r;
        forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2 - slater_term);
        if (factor_coul < static_cast<KK_FLOAT>(1.0))
          forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor*(static_cast<KK_FLOAT>(1.0)-slater_term);
        forcecoul *= r2inv;

        if (EVFLAG && eflag) {
          const KK_FLOAT eslater = (static_cast<KK_FLOAT>(1.0) + rlamda)*exp2rl;
          ecoul = prefactor*(erfc - eslater);
          if (factor_coul < static_cast<KK_FLOAT>(1.0))
            ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor*(static_cast<KK_FLOAT>(1.0)-eslater);
        }
      }

      const KK_FLOAT fpair = forcedpd + forcecoul;

      fx += fpair*delx;
      fy += fpair*dely;
      fz += fpair*delz;

      a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fpair*delx);
      a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fpair*dely);
      a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fpair*delz);

      if (EVFLAG) {
        if (eflag_global) {
          ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
          ev.ecoul += static_cast<KK_ACC_FLOAT>(ecoul);
        }
        if (eflag_atom || vflag_either)
          this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl+ecoul,fpair,delx,dely,delz);
      }
    }
  }
  a_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
  a_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
  a_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairDPDCoulSlaterLongKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
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
    const KK_FLOAT v0 = delx*delx*fpair;
    const KK_FLOAT v1 = dely*dely*fpair;
    const KK_FLOAT v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair;
    const KK_FLOAT v4 = delx*delz*fpair;
    const KK_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v0);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v1);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v2);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v3);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v4);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v5);
    }

    if (vflag_atom) {
      a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v0);
      a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v1);
      a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v2);
      a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v3);
      a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v4);
      a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v5);
      a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v0);
      a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v1);
      a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v2);
      a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v3);
      a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v4);
      a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v5);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDCoulSlaterLongKokkos<DeviceType>::allocate()
{
  PairDPDCoulSlaterLong::allocate();
  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType>("PairDPDCoulSlaterLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int PairDPDCoulSlaterLongKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDPDCoulSlaterLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairDPDCoulSlaterLong::init_one(i,j);

  k_params.view_host()(i,j).cut_dpd = static_cast<KK_FLOAT>(cut_dpd[i][j]);
  k_params.view_host()(i,j).cut_dpdsq = static_cast<KK_FLOAT>(cut_dpdsq[i][j]);
  k_params.view_host()(i,j).cut_slatersq = static_cast<KK_FLOAT>(cut_slatersq[i][j]);
  k_params.view_host()(i,j).a0 = static_cast<KK_FLOAT>(a0[i][j]);
  k_params.view_host()(i,j).gamma = static_cast<KK_FLOAT>(gamma[i][j]);
  k_params.view_host()(i,j).sigma = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  k_params.modify_host();

  k_cutsq.view_host()(i,j) = cutone*cutone;
  k_cutsq.view_host()(j,i) = k_cutsq.view_host()(i,j);
  k_cutsq.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairDPDCoulSlaterLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDPDCoulSlaterLongKokkos<LMPHostType>;
#endif
}
