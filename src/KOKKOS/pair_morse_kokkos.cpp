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
   Contributing authors: Stefan Paquay (Eindhoven University of Technology)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_morse_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMorseKokkos<DeviceType>::PairMorseKokkos(LAMMPS *lmp) : PairMorse(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMorseKokkos<DeviceType>::~PairMorseKokkos()
{
  if (allocated) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_ffloat_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;


  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  if (eflag || vflag) ev_setup(eflag,vflag,0);
  else evflag = vflag_fdotr = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memory->destroy_kokkos(k_vatom,vatom);
    memory->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairMorseKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairMorseKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT rr = sqrt(rsq);
  const F_FLOAT r0 = STACKPARAMS ? m_params[itype][jtype].r0 : params(itype,jtype).r0;
  const F_FLOAT d0 = STACKPARAMS ? m_params[itype][jtype].d0 : params(itype,jtype).d0;
  const F_FLOAT aa = STACKPARAMS ? m_params[itype][jtype].alpha : params(itype,jtype).alpha;
  const F_FLOAT dr = rr - r0;
  
  // U  =  d0 * [ exp( -2*a*(x-r0)) - 2*exp(-a*(x-r0)) ]
  // f  = -2*a*d0*[ -exp( -2*a*(x-r0) ) + exp( -a*(x-r0) ) ] * grad(r)
  //    = +2*a*d0*[  exp( -2*a*(x-r0) ) - exp( -a*(x-r0) ) ] * grad(r)
  const F_FLOAT dexp    = exp( -aa*dr );
  const F_FLOAT forcelj = 2*aa*d0*dexp*(dexp-1.0);

  return forcelj / rr;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairMorseKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT rr = sqrt(rsq);
  const F_FLOAT r0 = STACKPARAMS ? m_params[itype][jtype].r0 : params(itype,jtype).r0;
  const F_FLOAT d0 = STACKPARAMS ? m_params[itype][jtype].d0 : params(itype,jtype).d0;
  const F_FLOAT aa = STACKPARAMS ? m_params[itype][jtype].alpha : params(itype,jtype).alpha;
  const F_FLOAT dr = rr - r0;
  
  // U  =  d0 * [ exp( -2*a*(x-r0)) - 2*exp(-a*(x-r0)) ]
  // f  = -2*a*d0*[ -exp( -2*a*(x-r0) ) + exp( -a*(x-r0) ) ] * grad(r)
  //    = +2*a*d0*[  exp( -2*a*(x-r0) ) - exp( -a*(x-r0) ) ] * grad(r)
  const F_FLOAT dexp    = exp( -aa*dr );

  return d0 * dexp * ( dexp - 2.0 );
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::allocate()
{
  PairMorse::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memory->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_morse**,Kokkos::LayoutRight,DeviceType>("PairMorse::params",n+1,n+1);
  params = k_params.d_view;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairMorse::settings(1,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::init_style()
{
  PairMorse::init_style();

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
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else if (neighflag == N2) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with morse/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
// Rewrite this.
template<class DeviceType>
double PairMorseKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairMorse::init_one(i,j);

  k_params.h_view(i,j).d0     = d0[i][j];
  k_params.h_view(i,j).alpha  = alpha[i][j];
  k_params.h_view(i,j).r0     = r0[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cutsq  = cutone*cutone;
  k_params.h_view(j,i)        = k_params.h_view(i,j);
  
  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}



namespace LAMMPS_NS {
template class PairMorseKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairMorseKokkos<LMPHostType>;
#endif
}

