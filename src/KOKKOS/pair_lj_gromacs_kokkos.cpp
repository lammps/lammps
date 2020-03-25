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

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_gromacs_kokkos.h"
#include <cmath>
#include <cstring>
#include "kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "respa.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJGromacsKokkos<DeviceType>::PairLJGromacsKokkos(LAMMPS *lmp):PairLJGromacs(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
  cut_inner = NULL;
  cut_inner_sq = NULL;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJGromacsKokkos<DeviceType>::~PairLJGromacsKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_ffloat_2d();
    k_cut_inner_sq = DAT::tdual_ffloat_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
    cut_inner = NULL;
    cut_inner_sq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJGromacsKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  cut_inner = NULL;
  cut_inner_sq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJGromacsKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_inner.template sync<DeviceType>();
  k_cut_inner_sq.template sync<DeviceType>();
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
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev = pair_compute<PairLJGromacsKokkos<DeviceType>,CoulLongTable<0> >
      (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += 0.0;
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
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute LJ GROMACS pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJGromacsKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {

  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  if (rsq > (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT tlj = r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
    const F_FLOAT fswitch = r*tlj*tlj*
            ((STACKPARAMS?m_params[itype][jtype].ljsw1:params(itype,jtype).ljsw1) +
             (STACKPARAMS?m_params[itype][jtype].ljsw2:params(itype,jtype).ljsw2)*tlj);
    forcelj += fswitch;
  }
  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ GROMACS pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJGromacsKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {

  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT englj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4));
  englj += (STACKPARAMS?m_params[itype][jtype].ljsw5:params(itype,jtype).ljsw5);

  if (rsq > (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT tlj = r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
    const F_FLOAT eswitch = tlj*tlj*tlj *
            ((STACKPARAMS?m_params[itype][jtype].ljsw3:params(itype,jtype).ljsw3) +
             (STACKPARAMS?m_params[itype][jtype].ljsw4:params(itype,jtype).ljsw4)*tlj);
    englj += eswitch;
  }
  return englj;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJGromacsKokkos<DeviceType>::allocate()
{
  PairLJGromacs::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_inner);
  memoryKK->create_kokkos(k_cut_inner,cut_inner,n+1,n+1,"pair:cut_inner");
  d_cut_inner = k_cut_inner.template view<DeviceType>();

  memory->destroy(cut_inner_sq);
  memoryKK->create_kokkos(k_cut_inner_sq,cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  d_cut_inner_sq = k_cut_inner_sq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>("PairLJGromacs::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJGromacsKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  PairLJGromacs::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJGromacsKokkos<DeviceType>::init_style()
{
  PairLJGromacs::init_style();

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
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/gromacs/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJGromacsKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJGromacs::init_one(i,j);
  double cut_inner_sqm = cut_inner_sq[i][j];

  k_params.h_view(i,j).lj1 = lj1[i][j];
  k_params.h_view(i,j).lj2 = lj2[i][j];
  k_params.h_view(i,j).lj3 = lj3[i][j];
  k_params.h_view(i,j).lj4 = lj4[i][j];
  k_params.h_view(i,j).ljsw1 = ljsw1[i][j];
  k_params.h_view(i,j).ljsw2 = ljsw2[i][j];
  k_params.h_view(i,j).ljsw3 = ljsw3[i][j];
  k_params.h_view(i,j).ljsw4 = ljsw4[i][j];
  k_params.h_view(i,j).ljsw5 = ljsw5[i][j];
  k_params.h_view(i,j).cut_inner_sq = cut_inner_sqm;
  k_params.h_view(i,j).cut_inner = cut_inner[i][j];

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_inner_sq[j][i] = m_cut_inner_sq[i][j] = cut_inner_sqm;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_cut_inner_sq.h_view(i,j) = k_cut_inner_sq.h_view(j,i) = cut_inner_sqm;
  k_cut_inner_sq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJGromacsKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairLJGromacsKokkos<LMPHostType>;
#endif
}

