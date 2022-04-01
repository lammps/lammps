// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#include "pair_dpd_ext_kokkos.h"

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

#define EPSILON 1.0e-10


template<class DeviceType>
PairDPDExtKokkos<DeviceType>::PairDPDExtKokkos(class LAMMPS *lmp) :
  PairDPDExt(lmp) ,
#ifdef DPD_USE_RAN_MARS
  rand_pool(0 /* unused */, lmp)
#else
  rand_pool()
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDExtKokkos<DeviceType>::~PairDPDExtKokkos() {
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
void PairDPDExtKokkos<DeviceType>::init_style()
{
  PairDPDExt::init_style();

#ifdef DPD_USE_RAN_MARS
  rand_pool.init(random,seed);
#else
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  rand_pool.init(seed + comm->me,unique_token.size());
#endif

  neighflag = lmp->kokkos->neighflag;

  if (force->newton_pair == 0 || neighflag == FULL )
    error->all(FLERR,"Must use half neighbor list style and newton on with pair dpd/ext/kk");

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);

  if (neighflag == FULL)
    request->enable_full();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDExtKokkos<DeviceType>::compute(int eflagin, int vflagin)
{
  eflag = eflagin; vflag = vflagin;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag);

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

  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  dtinvsqrt = 1.0/sqrt(update->dt);

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // loop over neighbors of my atoms

  int inum = list->inum;
  EV_FLOAT ev;
  copymode = 1;
  if (neighflag == HALF) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALF,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALF,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALF,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALF,0,0> >(0,inum),*this);
    }
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALFTHREAD,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALFTHREAD,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALFTHREAD,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<HALFTHREAD,0,0> >(0,inum),*this);
    }
  } else if (neighflag == FULL) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<FULL,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<FULL,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<FULL,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDPDExtKokkos<FULL,0,0> >(0,inum),*this);
    }
  }

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
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;

  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  else atomKK->modified(execution_space,F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairDPDExtKokkos<DeviceType>::operator() (TagDPDExtKokkos<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagDPDExtKokkos<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairDPDExtKokkos<DeviceType>::operator() (TagDPDExtKokkos<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT &ev) const {

  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;

  int i,j,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpairx,fpairy,fpairz,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wdPar,wdPerp,randnum,randnumx,randnumy,randnumz,factor_dpd;
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
    double P[3][3];
    j = d_neighbors(i,jj);
    factor_dpd = special_lj[sbmask(j)];
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

      P[0][0] = 1.0 - delx*delx*rinv*rinv;
      P[0][1] =     - delx*dely*rinv*rinv;
      P[0][2] =     - delx*delz*rinv*rinv;

      P[1][0] = P[0][1];
      P[1][1] = 1.0 - dely*dely*rinv*rinv;
      P[1][2] =     - dely*delz*rinv*rinv;

      P[2][0] = P[0][2];
      P[2][1] = P[1][2];
      P[2][2] = 1.0 - delz*delz*rinv*rinv;

      wd = 1.0 - r/params(itype,jtype).cut;
      wdPar = pow(wd,params(itype,jtype).ws);
      wdPerp = pow(wd,params(itype,jtype).wsT);

      randnum  = rand_gen.normal();
      randnumx = rand_gen.normal();
      randnumy = rand_gen.normal();
      randnumz = rand_gen.normal();
      // conservative force
      fpair = params(itype,jtype).a0*wd;

      // drag force - parallel
      fpair -= params(itype,jtype).gamma*wdPar*wdPar*dot*rinv;

      // random force - parallel
      fpair += params(itype,jtype).sigma*wdPar*randnum*dtinvsqrt;

      fpairx = fpair*rinv*delx;
      fpairy = fpair*rinv*dely;
      fpairz = fpair*rinv*delz;

      // drag force - perpendicular
      fpairx -= params(itype,jtype).gammaT*wdPerp*wdPerp*
                (P[0][0]*delvx + P[0][1]*delvy + P[0][2]*delvz);
      fpairy -= params(itype,jtype).gammaT*wdPerp*wdPerp*
                (P[1][0]*delvx + P[1][1]*delvy + P[1][2]*delvz);
      fpairz -= params(itype,jtype).gammaT*wdPerp*wdPerp*
                (P[2][0]*delvx + P[2][1]*delvy + P[2][2]*delvz);

      // random force - perpendicular
      fpairx += params(itype,jtype).sigmaT*wdPerp*
                (P[0][0]*randnumx + P[0][1]*randnumy + P[0][2]*randnumz)*dtinvsqrt;
      fpairy += params(itype,jtype).sigmaT*wdPerp*
                (P[1][0]*randnumx + P[1][1]*randnumy + P[1][2]*randnumz)*dtinvsqrt;
      fpairz += params(itype,jtype).sigmaT*wdPerp*
                (P[2][0]*randnumx + P[2][1]*randnumy + P[2][2]*randnumz)*dtinvsqrt;

      fpairx *= factor_dpd;
      fpairy *= factor_dpd;
      fpairz *= factor_dpd;

      fx += fpairx;
      fy += fpairy;
      fz += fpairz;
      if ((neighflag==HALF || neighflag==HALFTHREAD) && (NEWTON_PAIR || j < nlocal) ) {
        a_f(j,0) -= fpairx;
        a_f(j,1) -= fpairy;
        a_f(j,2) -= fpairz;
      }

      if (eflag) {
        // unshifted eng of conservative term:
        // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
        // eng shifted to 0.0 at cutoff
        evdwl = 0.5*params(itype,jtype).a0*params(itype,jtype).cut* wd*wd;
        evdwl *= factor_dpd;
        if (EVFLAG)
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      }
      if (EVFLAG)
	this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }
  a_f(i,0) += fx;
  a_f(i,1) += fy;
  a_f(i,2) += fz;
  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairDPDExtKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_vatom = k_vatom.view<DeviceType>();

  if (EFLAG) {
    if (eflag_atom) {
      const E_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
        if (NEWTON_PAIR || j < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
      } else {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          a_vatom(i,0) += 0.5*v0;
          a_vatom(i,1) += 0.5*v1;
          a_vatom(i,2) += 0.5*v2;
          a_vatom(i,3) += 0.5*v3;
          a_vatom(i,4) += 0.5*v4;
          a_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        a_vatom(j,0) += 0.5*v0;
        a_vatom(j,1) += 0.5*v1;
        a_vatom(j,2) += 0.5*v2;
        a_vatom(j,3) += 0.5*v3;
        a_vatom(j,4) += 0.5*v4;
        a_vatom(j,5) += 0.5*v5;
        }
      } else {
        a_vatom(i,0) += 0.5*v0;
        a_vatom(i,1) += 0.5*v1;
        a_vatom(i,2) += 0.5*v2;
        a_vatom(i,3) += 0.5*v3;
        a_vatom(i,4) += 0.5*v4;
        a_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDExtKokkos<DeviceType>::allocate()
{
  PairDPDExt::allocate();
  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType>("PairDPDExt::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairDPDExtKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDPDExtKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairDPDExt::init_one(i,j);

  k_params.h_view(i,j).cut = cut[i][j];
  k_params.h_view(i,j).ws = ws[i][j];
  k_params.h_view(i,j).wsT = wsT[i][j];
  k_params.h_view(i,j).a0 = a0[i][j];
  k_params.h_view(i,j).gamma = gamma[i][j];
  k_params.h_view(i,j).sigma = sigma[i][j];
  k_params.h_view(i,j).gammaT = gammaT[i][j];
  k_params.h_view(i,j).sigmaT = sigmaT[i][j];
  k_params.h_view(j,i) = k_params.h_view(i,j);

  k_params.template modify<LMPHostType>();

  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairDPDExtKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDPDExtKokkos<LMPHostType>;
#endif
}
