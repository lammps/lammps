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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_buck_kokkos.h"
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
PairBuckKokkos<DeviceType>::PairBuckKokkos(LAMMPS *lmp) : PairBuck(lmp)
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
PairBuckKokkos<DeviceType>::~PairBuckKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_ffloat_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBuckKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  copymode = 1;
  EV_FLOAT ev = pair_compute<PairBuckKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairBuckKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  const F_FLOAT r = sqrt(rsq);
  const F_FLOAT rexp = exp(-r*(STACKPARAMS?m_params[itype][jtype].rhoinv:params(itype,jtype).rhoinv));

  const F_FLOAT forcebuck =
     (STACKPARAMS?m_params[itype][jtype].buck1:params(itype,jtype).buck1)*r*rexp -
     (STACKPARAMS?m_params[itype][jtype].buck2:params(itype,jtype).buck2)*r6inv;

  return forcebuck*r2inv;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairBuckKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  const F_FLOAT r = sqrt(rsq);
  const F_FLOAT rexp = exp(-r*(STACKPARAMS?m_params[itype][jtype].rhoinv:params(itype,jtype).rhoinv));


  return (STACKPARAMS?m_params[itype][jtype].a:params(itype,jtype).a)*rexp -
                (STACKPARAMS?m_params[itype][jtype].c:params(itype,jtype).c)*r6inv -
                (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuckKokkos<DeviceType>::allocate()
{
  PairBuck::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_buck**,Kokkos::LayoutRight,DeviceType>("PairBuck::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuckKokkos<DeviceType>::init_style()
{
  PairBuck::init_style();

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
    error->all(FLERR,"Cannot use chosen neighbor list style with buck/kk");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairBuckKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBuck::init_one(i,j);

  k_params.h_view(i,j).a = a[i][j];
  k_params.h_view(i,j).c = c[i][j];
  k_params.h_view(i,j).rhoinv = rhoinv[i][j];
  k_params.h_view(i,j).buck1 = buck1[i][j];
  k_params.h_view(i,j).buck2 = buck2[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cutsq = cutone*cutone;
  k_params.h_view(j,i) = k_params.h_view(i,j);
  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }
  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairBuckKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairBuckKokkos<LMPHostType>;
#endif
}

