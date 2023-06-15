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
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "pair_dpd_fdt_energy_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDfdtEnergyKokkos<DeviceType>::PairDPDfdtEnergyKokkos(LAMMPS *lmp) :
  PairDPDfdtEnergy(lmp),
#ifdef DPD_USE_RAN_MARS
  rand_pool(0 /* unused */, lmp)
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
PairDPDfdtEnergyKokkos<DeviceType>::~PairDPDfdtEnergyKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  if (allocated) {
    memoryKK->destroy_kokkos(k_duCond,duCond);
    memoryKK->destroy_kokkos(k_duMech,duMech);
  }

  memoryKK->destroy_kokkos(k_cutsq,cutsq);

#ifdef DPD_USE_RAN_MARS
  rand_pool.destroy();
#endif
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::init_style()
{
  PairDPDfdtEnergy::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();

#ifdef DPD_USE_RAN_MARS
  rand_pool.init(random,seed);
#else
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  rand_pool.init(seed + comm->me,unique_token.size());
#endif
}

#if (defined(KOKKOS_ENABLE_CUDA) && defined(__CUDACC__)) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
// CUDA specialization of init_style to properly call rand_pool.init()
template<>
void PairDPDfdtEnergyKokkos<LMPDeviceSpace>::init_style()
{
  PairDPDfdtEnergy::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(0);
  request->set_kokkos_device(1);
  if (neighflag == FULL) request->enable_full();

#ifdef DPD_USE_RAN_MARS
  rand_pool.init(random,seed);
#else
  rand_pool.init(seed + comm->me,4*32768 /*fake max_hardware_threads()*/);
#endif
}
#endif

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

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

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  dpdTheta = atomKK->k_dpdTheta.view<DeviceType>();

  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  atomKK->sync(execution_space,X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK);
  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  else atomKK->modified(execution_space,F_MASK);

  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;
  dtinvsqrt = 1.0/sqrt(update->dt);

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  boltz = force->boltz;
  ftm2v = force->ftm2v;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (splitFDT_flag) {
    if (!a0_is_zero) {
      if (atom->ntypes > MAX_TYPES_STACKPARAMS) {
        if (neighflag == HALF) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,1,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,1,0,false> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,0,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,0,0,false> >(0,inum),*this);
          }
        } else if (neighflag == HALFTHREAD) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,1,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,1,0,false> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,0,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,0,0,false> >(0,inum),*this);
          }
        } else if (neighflag == FULL) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,1,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,1,0,false> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,0,1,false> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,0,0,false> >(0,inum),*this);
          }
        }
      } else {
        if (neighflag == HALF) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,1,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,1,0,true> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,0,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALF,0,0,true> >(0,inum),*this);
          }
        } else if (neighflag == HALFTHREAD) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,1,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,1,0,true> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,0,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<HALFTHREAD,0,0,true> >(0,inum),*this);
          }
        } else if (neighflag == FULL) {
          if (newton_pair) {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,1,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,1,0,true> >(0,inum),*this);
          } else {
            if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,0,1,true> >(0,inum),*this,ev);
            else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeSplit<FULL,0,0,true> >(0,inum),*this);
          }
        }
      }
    }
  } else {

    // Allocate memory for duCond and duMech
    if (allocated) {
      memoryKK->destroy_kokkos(k_duCond,duCond);
      memoryKK->destroy_kokkos(k_duMech,duMech);
    }
    memoryKK->create_kokkos(k_duCond,duCond,nlocal+nghost,"pair:duCond");
    memoryKK->create_kokkos(k_duMech,duMech,nlocal+nghost,"pair:duMech");
    d_duCond = k_duCond.view<DeviceType>();
    d_duMech = k_duMech.view<DeviceType>();
    h_duCond = k_duCond.h_view;
    h_duMech = k_duMech.h_view;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyZero>(0,nlocal+nghost),*this);

    atomKK->sync(execution_space,V_MASK | DPDTHETA_MASK | RMASS_MASK);
    atomKK->k_mass.sync<DeviceType>();

    // loop over neighbors of my atoms

    if (atom->ntypes > MAX_TYPES_STACKPARAMS) {
      if (neighflag == HALF) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,0,0,false> >(0,inum),*this);
        }
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,0,0,false> >(0,inum),*this);
        }
      } else if (neighflag == FULL) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,0,0,false> >(0,inum),*this);
        }
      }
    } else {
      if (neighflag == HALF) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALF,0,0,false> >(0,inum),*this);
        }
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<HALFTHREAD,0,0,false> >(0,inum),*this);
        }
      } else if (neighflag == FULL) {
        if (newton_pair) {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,1,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,1,0,false> >(0,inum),*this);
        } else {
          if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,0,1,false> >(0,inum),*this,ev);
          else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairDPDfdtEnergyComputeNoSplit<FULL,0,0,false> >(0,inum),*this);
        }
      }
    }

    // Communicate the ghost delta energies to the locally owned atoms

    // this memory transfer can be removed when fix_dpd_fdt_energy_kokkos is added
    k_duCond.template modify<DeviceType>();
    k_duCond.template sync<LMPHostType>();
    k_duMech.template modify<DeviceType>();
    k_duMech.template sync<LMPHostType>();
    comm->reverse_comm(this);
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
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::operator()(TagPairDPDfdtEnergyZero, const int &ii) const {
  d_duCond[ii] = 0.0;
  d_duMech[ii] = 0.0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int &ii, EV_FLOAT& ev) const {

  // The f array is atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;

  int i,j,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rinv,wd,wr,factor_dpd;

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;

  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors(i,jj);
    factor_dpd = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    double cutsq_ij = STACKPARAMS?m_cutsq[itype][jtype]:d_cutsq(itype,jtype);
    if (rsq < cutsq_ij) {
      r = sqrt(rsq);
      if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
      rinv = 1.0/r;
      double cut_ij = STACKPARAMS?m_params[itype][jtype].cut:params(itype,jtype).cut;
      wr = 1.0 - r/cut_ij;
      wd = wr*wr;

      // conservative force = a0 * wr
      double a0_ij = STACKPARAMS?m_params[itype][jtype].a0:params(itype,jtype).a0;
      fpair = a0_ij*wr;
      fpair *= factor_dpd*rinv;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= delx*fpair;
        a_f(j,1) -= dely*fpair;
        a_f(j,2) -= delz*fpair;
      }

      if (eflag) {
        // unshifted eng of conservative term:
        // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/d_cut(itype,jtype));
        // eng shifted to 0.0 at cutoff
        evdwl = 0.5*a0_ij*cut_ij * wd;
        evdwl *= factor_dpd;
        if (EVFLAG)
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      }

      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int &ii, EV_FLOAT& ev) const {

  // These array are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_duCond = d_duCond;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_duMech = d_duMech;

  int i,j,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,wd,wr,factor_dpd,uTmp;
  double dot,randnum;

  double kappa_ij, alpha_ij, theta_ij, gamma_ij;
  double mass_i, mass_j;
  double massinv_i, massinv_j;
  double randPair, mu_ij;

  rand_type rand_gen = rand_pool.get_state();

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  vxtmp = v(i,0);
  vytmp = v(i,1);
  vztmp = v(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;

  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors(i,jj);
    factor_dpd = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    double cutsq_ij = STACKPARAMS?m_cutsq[itype][jtype]:d_cutsq(itype,jtype);
    if (rsq < cutsq_ij) {
      r = sqrt(rsq);
      if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
      rinv = 1.0/r;
      double cut_ij = STACKPARAMS?m_params[itype][jtype].cut:params(itype,jtype).cut;
      wr = 1.0 - r/cut_ij;
      wd = wr*wr;

      delvx = vxtmp - v(j,0);
      delvy = vytmp - v(j,1);
      delvz = vztmp - v(j,2);
      dot = delx*delvx + dely*delvy + delz*delvz;
      randnum = rand_gen.normal();

      // Compute the current temperature
      theta_ij = 0.5*(1.0/dpdTheta[i] + 1.0/dpdTheta[j]);
      theta_ij = 1.0/theta_ij;

      double sigma_ij = STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma;
      gamma_ij = sigma_ij*sigma_ij
                 / (2.0*boltz*theta_ij);

      // conservative force = a0 * wr
      // drag force = -gamma * wr^2 * (delx dot delv) / r
      // random force = sigma * wr * rnd * dtinvsqrt;

      double a0_ij = STACKPARAMS?m_params[itype][jtype].a0:params(itype,jtype).a0;
      fpair = a0_ij*wr;
      fpair -= gamma_ij*wd*dot*rinv;
      fpair += sigma_ij*wr*randnum*dtinvsqrt;
      fpair *= factor_dpd*rinv;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= delx*fpair;
        a_f(j,1) -= dely*fpair;
        a_f(j,2) -= delz*fpair;
      }

      if (rmass.data()) {
        mass_i = rmass[i];
        mass_j = rmass[j];
      } else {
        mass_i = mass[itype];
        mass_j = mass[jtype];
      }
      massinv_i = 1.0 / mass_i;
      massinv_j = 1.0 / mass_j;

      // Compute the mechanical and conductive energy, uMech and uCond
      mu_ij = massinv_i + massinv_j;
      mu_ij *= ftm2v;

      uTmp = gamma_ij*wd*rinv*rinv*dot*dot
             - 0.5*sigma_ij*sigma_ij*mu_ij*wd;
      uTmp -= sigma_ij*wr*rinv*dot*randnum*dtinvsqrt;
      uTmp *= 0.5;

      a_duMech[i] += uTmp;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_duMech[j] += uTmp;
      }

      // Compute uCond
      randnum = rand_gen.normal();
      kappa_ij = STACKPARAMS?m_params[itype][jtype].kappa:params(itype,jtype).kappa;
      alpha_ij = STACKPARAMS?m_params[itype][jtype].alpha:params(itype,jtype).alpha;
      randPair = alpha_ij*wr*randnum*dtinvsqrt;

      uTmp = kappa_ij*(1.0/dpdTheta[i] - 1.0/dpdTheta[j])*wd;
      uTmp += randPair;

      a_duCond[i] += uTmp;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_duCond[j] -= uTmp;
      }

      if (eflag) {
        // unshifted eng of conservative term:
        // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/d_cut(itype,jtype));
        // eng shifted to 0.0 at cutoff
        evdwl = 0.5*a0_ij*cut_ij * wd;
        evdwl *= factor_dpd;
        if (EVFLAG)
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      }

      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;

  rand_pool.free_state(rand_gen);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>(), ii, ev);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::allocate()
{
  PairDPDfdtEnergy::allocate();

  int n = atom->ntypes;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType>("PairDPDfdtEnergy::params",n+1,n+1);
  params = k_params.template view<DeviceType>();

  if (!splitFDT_flag) {
    memory->destroy(duCond);
    memory->destroy(duMech);
    memoryKK->create_kokkos(k_duCond,duCond,nlocal+nghost+1,"pair:duCond");
    memoryKK->create_kokkos(k_duMech,duMech,nlocal+nghost+1,"pair:duMech");
    d_duCond = k_duCond.view<DeviceType>();
    d_duMech = k_duMech.view<DeviceType>();
    h_duCond = k_duCond.h_view;
    h_duMech = k_duMech.h_view;
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDPDfdtEnergyKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairDPDfdtEnergy::init_one(i,j);

  k_params.h_view(i,j).cut = cut[i][j];
  k_params.h_view(i,j).a0 = a0[i][j];
  k_params.h_view(i,j).sigma = sigma[i][j];
  k_params.h_view(i,j).kappa = kappa[i][j];
  k_params.h_view(i,j).alpha = alpha[i][j];
  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairDPDfdtEnergyKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  if (EFLAG) {
    if (eflag_atom) {
      const E_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) v_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) v_eatom[j] += epairhalf;
      } else {
        v_eatom[i] += epairhalf;
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
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
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
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
        if (NEWTON_PAIR || j < nlocal) {
          v_vatom(j,0) += 0.5*v0;
          v_vatom(j,1) += 0.5*v1;
          v_vatom(j,2) += 0.5*v2;
          v_vatom(j,3) += 0.5*v3;
          v_vatom(j,4) += 0.5*v4;
          v_vatom(j,5) += 0.5*v5;
        }
      } else {
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairDPDfdtEnergyKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairDPDfdtEnergyKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDPDfdtEnergyKokkos<LMPHostType>;
#endif
}
