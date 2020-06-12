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

#include "pair_coul_cut_kokkos.h"
#include <cmath>
#include "kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairCoulCutKokkos<Space>::PairCoulCutKokkos(LAMMPS *lmp) : PairCoulCut(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairCoulCutKokkos<Space>::~PairCoulCutKokkos()
{
  if (allocated)
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulCutKokkos<Space>::cleanup_copy() {
 // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulCutKokkos<Space>::compute(int eflag_in, int vflag_in)
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
  DualViewHelper<Space>::sync(k_cut_ljsq);
  DualViewHelper<Space>::sync(k_cut_coulsq);
  DualViewHelper<Space>::sync(k_params);
  if (eflag || vflag) atomKK->modified(Space,datamask_modify);
  else atomKK->modified(Space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  q = DualViewHelper<Space>::view(atomKK->k_q);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<Space,PairCoulCutKokkos<Space>,void >
    (this,(NeighListKokkos<Space>*)list);

  if (eflag) eng_coul += ev.ecoul;
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
KK_FLOAT PairCoulCutKokkos<Space>::
compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT r2inv = 1.0/rsq;
  const KK_FLOAT rinv = sqrt(r2inv);
  KK_FLOAT forcecoul;

  forcecoul = qqrd2e*(STACKPARAMS?m_params[itype][jtype].scale:params(itype,jtype).scale)*
    qtmp *q(j) *rinv;

  return factor_coul*forcecoul*r2inv;
}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulCutKokkos<Space>::
compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT r2inv = 1.0/rsq;
  const KK_FLOAT rinv = sqrt(r2inv);

  return factor_coul*qqrd2e * (STACKPARAMS?m_params[itype][jtype].scale:params(itype,jtype).scale)
    * qtmp *q(j)*rinv;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulCutKokkos<Space>::allocate()
{
  PairCoulCut::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = DualViewHelper<Space>::view(k_cutsq);

  k_cut_ljsq = DAT::tdual_float_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_ljsq = DualViewHelper<Space>::view(k_cut_ljsq);
  k_cut_coulsq = DAT::tdual_float_2d("pair:cut_coulsq",n+1,n+1);
  d_cut_coulsq = DualViewHelper<Space>::view(k_cut_coulsq);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulCut::params",n+1,n+1);
  params = DualViewHelper<Space>::view(k_params);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulCutKokkos<Space>::settings(int narg, char **arg)
{
  // \todo check what should be the limit on narg
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairCoulCut::settings(1,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairCoulCutKokkos<Space>::init_style()
{
  PairCoulCut::init_style();

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
    error->all(FLERR,"Cannot use chosen neighbor list style with coul/cut/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairCoulCutKokkos<Space>::init_one(int i, int j)
{
  KK_FLOAT cutone = PairCoulCut::init_one(i,j);

  k_params.h_view(i,j).scale = scale[i][j];
  k_params.h_view(i,j).cutsq = cutone*cutone;
  k_params.h_view(j,i) = k_params.h_view(i,j);

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cutone*cutone;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cutone*cutone;
  }
  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.h_view(i,j) = cutone*cutone;
  k_cut_ljsq.modify_host();
  k_cut_coulsq.h_view(i,j) = cutone*cutone;
  k_cut_coulsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairCoulCutKokkos<Device>;
template class PairCoulCutKokkos<Host>;
}

