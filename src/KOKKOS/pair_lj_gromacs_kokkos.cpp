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

template<ExecutionSpace Space>
PairLJGromacsKokkos<Space>::PairLJGromacsKokkos(LAMMPS *lmp):PairLJGromacs(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
  cut_inner = NULL;
  cut_inner_sq = NULL;

}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairLJGromacsKokkos<Space>::~PairLJGromacsKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_float_2d();
    k_cut_inner_sq = DAT::tdual_float_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
    cut_inner = NULL;
    cut_inner_sq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJGromacsKokkos<Space>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  cut_inner = NULL;
  cut_inner_sq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJGromacsKokkos<Space>::compute(int eflag_in, int vflag_in)
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

  atomKK->sync(Space,datamask_read);
  DualViewHelper<Space>::sync(k_cutsq);
  DualViewHelper<Space>::sync(k_cut_inner);
  DualViewHelper<Space>::sync(k_cut_inner_sq);
  DualViewHelper<Space>::sync(k_params);
  if (eflag || vflag) atomKK->modified(Space,datamask_modify);
  else atomKK->modified(Space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  q = DualViewHelper<Space>::view(atomKK->k_q);
  type = DualViewHelper<Space>::view(atomKK->k_type);
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

  EV_FLOAT ev = pair_compute<Space,PairLJGromacsKokkos<Space>,CoulLongTable<0> >
      (this,(NeighListKokkos<Space>*)list);

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
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute LJ GROMACS pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJGromacsKokkos<Space>::
compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {

  const KK_FLOAT r2inv = 1.0/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  KK_FLOAT forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  if (rsq > (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT tlj = r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
    const KK_FLOAT fswitch = r*tlj*tlj*
            ((STACKPARAMS?m_params[itype][jtype].ljsw1:params(itype,jtype).ljsw1) +
             (STACKPARAMS?m_params[itype][jtype].ljsw2:params(itype,jtype).ljsw2)*tlj);
    forcelj += fswitch;
  }
  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ GROMACS pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJGromacsKokkos<Space>::
compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {

  const KK_FLOAT r2inv = 1.0/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  KK_FLOAT englj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4));
  englj += (STACKPARAMS?m_params[itype][jtype].ljsw5:params(itype,jtype).ljsw5);

  if (rsq > (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT tlj = r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
    const KK_FLOAT eswitch = tlj*tlj*tlj *
            ((STACKPARAMS?m_params[itype][jtype].ljsw3:params(itype,jtype).ljsw3) +
             (STACKPARAMS?m_params[itype][jtype].ljsw4:params(itype,jtype).ljsw4)*tlj);
    englj += eswitch;
  }
  return englj;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJGromacsKokkos<Space>::allocate()
{
  PairLJGromacs::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = DualViewHelper<Space>::view(k_cutsq);

  memory->destroy(cut_inner);
  memoryKK->create_kokkos(k_cut_inner,cut_inner,n+1,n+1,"pair:cut_inner");
  d_cut_inner = DualViewHelper<Space>::view(k_cut_inner);

  memory->destroy(cut_inner_sq);
  memoryKK->create_kokkos(k_cut_inner_sq,cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  d_cut_inner_sq = DualViewHelper<Space>::view(k_cut_inner_sq);

  k_params = Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>("PairLJGromacs::params",n+1,n+1);
  params = DualViewHelper<Space>::view(k_params);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJGromacsKokkos<Space>::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  PairLJGromacs::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairLJGromacsKokkos<Space>::init_style()
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
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/gromacs/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairLJGromacsKokkos<Space>::init_one(int i, int j)
{
  KK_FLOAT cutone = PairLJGromacs::init_one(i,j);
  KK_FLOAT cut_inner_sqm = cut_inner_sq[i][j];

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
  k_cutsq.modify_host();
  k_cut_inner_sq.h_view(i,j) = k_cut_inner_sq.h_view(j,i) = cut_inner_sqm;
  k_cut_inner_sq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJGromacsKokkos<Device>;
template class PairLJGromacsKokkos<Host>;
}

