/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_sdk_kokkos.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include "kokkos.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "respa.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

#include "lj_sdk_common.h"

using namespace LAMMPS_NS;
using namespace LJSDKParms;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairLJSDKKokkos<Space>::PairLJSDKKokkos(LAMMPS *lmp) : PairLJSDK(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairLJSDKKokkos<Space>::~PairLJSDKKokkos()
{
  if (allocated) {
    k_cutsq = DAT::tdual_float_2d();
    memory->sfree(cutsq);
    cutsq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJSDKKokkos<Space>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJSDKKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;


  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = DualViewHelper<Space>::view(k_eatom);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = DualViewHelper<Space>::view(k_vatom);
  }

  atomKK->sync(execution_space,datamask_read);
  DualViewHelper<Space>::sync(k_cutsq);
  DualViewHelper<Space>::sync(k_params);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  tag = DualViewHelper<Space>::view(atomKK->k_tag);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<Space,PairLJSDKKokkos<Space>,void >(this,(NeighListKokkos<Space>*)list);

  if (eflag) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSDKKokkos<Space>::
compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const KK_FLOAT r2inv = 1.0/rsq;
  const int ljt = (STACKPARAMS?m_params[itype][jtype].lj_type:params(itype,jtype).lj_type);

  const KK_FLOAT lj_1 =  (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj_2 =  (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2);

  /*if (ljt == LJ12_4) {

    const KK_FLOAT r4inv=r2inv*r2inv;
    return r4inv*(lj_1*r4inv*r4inv - lj_2) * r2inv;

  } else if (ljt == LJ9_6) {

    const KK_FLOAT r3inv = r2inv*sqrt(r2inv);
    const KK_FLOAT r6inv = r3inv*r3inv;
    return r6inv*(lj_1*r3inv - lj_2) * r2inv;

  } else if (ljt == LJ12_6) {

    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    return r6inv*(lj_1*r6inv - lj_2) * r2inv;

  }
  if(ljt!=LJ12_4 && ljt!=LJ9_6 && ljt!=LJ12_6) return 0.0;*/
  const KK_FLOAT r4inv=r2inv*r2inv;
  const KK_FLOAT r6inv=r2inv*r4inv;
  const KK_FLOAT a = ljt==LJ12_4?r4inv:r6inv;
  const KK_FLOAT b = ljt==LJ12_4?r4inv:(ljt==LJ9_6?1.0/sqrt(r2inv):r2inv);
  return a* ( lj_1*r6inv*b - lj_2 * r2inv);
}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSDKKokkos<Space>::
compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const KK_FLOAT r2inv = 1.0/rsq;
  const int ljt = (STACKPARAMS?m_params[itype][jtype].lj_type:params(itype,jtype).lj_type);

  const KK_FLOAT lj_3 =  (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3);
  const KK_FLOAT lj_4 =  (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4);
  const KK_FLOAT offset =  (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);

  if (ljt == LJ12_4) {
    const KK_FLOAT r4inv=r2inv*r2inv;

    return r4inv*(lj_3*r4inv*r4inv - lj_4) - offset;

  } else if (ljt == LJ9_6) {
    const KK_FLOAT r3inv = r2inv*sqrt(r2inv);
    const KK_FLOAT r6inv = r3inv*r3inv;
    return r6inv*(lj_3*r3inv - lj_4) - offset;

  } else if (ljt == LJ12_6) {
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    return r6inv*(lj_3*r6inv - lj_4) - offset;
  } else
    return 0.0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJSDKKokkos<Space>::allocate()
{
  PairLJSDK::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = DualViewHelper<Space>::view(k_cutsq);
  k_params = Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>("PairLJSDK::params",n+1,n+1);
  params = DualViewHelper<Space>::view(k_params);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJSDKKokkos<Space>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairLJSDK::settings(1,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJSDKKokkos<Space>::init_style()
{
  PairLJSDK::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/sdk/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairLJSDKKokkos<Space>::init_one(int i, int j)
{
  KK_FLOAT cutone = PairLJSDK::init_one(i,j);

  k_params.h_view(i,j).lj1 = lj1[i][j];
  k_params.h_view(i,j).lj2 = lj2[i][j];
  k_params.h_view(i,j).lj3 = lj3[i][j];
  k_params.h_view(i,j).lj4 = lj4[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cutsq = cutone*cutone;
  k_params.h_view(i,j).lj_type = lj_type[i][j];
  k_params.h_view(j,i) = k_params.h_view(i,j);
  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairLJSDKKokkos<Device>;
template class PairLJSDKKokkos<Host>;
}

