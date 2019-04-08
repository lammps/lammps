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
   Contributing authors: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_zbl_kokkos.h"
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
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "kokkos.h"

// From J.F. Zeigler, J. P. Biersack and U. Littmark,
// "The Stopping and Range of Ions in Matter" volume 1, Pergamon, 1985.

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace PairZBLConstants;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairZBLKokkos<DeviceType>::PairZBLKokkos(LAMMPS *lmp) : PairZBL(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairZBLKokkos<DeviceType>::~PairZBLKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
  }
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairZBLKokkos<DeviceType>::init_style()
{
  PairZBL::init_style();

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
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }

  Kokkos::deep_copy(d_cutsq,cut_globalsq);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairZBLKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  k_z.sync<DeviceType>();
  k_d1a.sync<DeviceType>();
  k_d2a.sync<DeviceType>();
  k_d3a.sync<DeviceType>();
  k_d4a.sync<DeviceType>();
  k_zze.sync<DeviceType>();
  k_sw1.sync<DeviceType>();
  k_sw2.sync<DeviceType>();
  k_sw3.sync<DeviceType>();
  k_sw4.sync<DeviceType>();
  k_sw5.sync<DeviceType>();

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairZBLKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += ev.evdwl;
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
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairZBLKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r = sqrt(rsq);
  F_FLOAT fpair = dzbldr(r, itype, jtype);

  if (rsq > cut_innersq) {
    const F_FLOAT t = r - cut_inner;
    const F_FLOAT fswitch = t*t *
           (d_sw1(itype,jtype) + d_sw2(itype,jtype)*t);
    fpair += fswitch;
  }

  fpair *= -1.0/r;
  return fpair;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairZBLKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r = sqrt(rsq);
  F_FLOAT evdwl = e_zbl(r, itype, jtype);
  evdwl += d_sw5(itype,jtype);
  if (rsq > cut_innersq) {
    const F_FLOAT t = r - cut_inner;
    const F_FLOAT eswitch = t*t*t *
      (d_sw3(itype,jtype) + d_sw4(itype,jtype)*t);
    evdwl += eswitch;
  }
  return evdwl;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairZBLKokkos<DeviceType>::allocate()
{
  PairZBL::allocate();

  int n = atom->ntypes;

  k_z   = DAT::tdual_ffloat_1d("pair_zbl:z  ",n+1);
  k_d1a = DAT::tdual_ffloat_2d_dl("pair_zbl:d1a",n+1,n+1);
  k_d2a = DAT::tdual_ffloat_2d_dl("pair_zbl:d2a",n+1,n+1);
  k_d3a = DAT::tdual_ffloat_2d_dl("pair_zbl:d3a",n+1,n+1);
  k_d4a = DAT::tdual_ffloat_2d_dl("pair_zbl:d4a",n+1,n+1);
  k_zze = DAT::tdual_ffloat_2d_dl("pair_zbl:zze",n+1,n+1);
  k_sw1 = DAT::tdual_ffloat_2d_dl("pair_zbl:sw1",n+1,n+1);
  k_sw2 = DAT::tdual_ffloat_2d_dl("pair_zbl:sw2",n+1,n+1);
  k_sw3 = DAT::tdual_ffloat_2d_dl("pair_zbl:sw3",n+1,n+1);
  k_sw4 = DAT::tdual_ffloat_2d_dl("pair_zbl:sw4",n+1,n+1);
  k_sw5 = DAT::tdual_ffloat_2d_dl("pair_zbl:sw5",n+1,n+1);

  d_z   = k_z.view<DeviceType>();
  d_d1a = k_d1a.view<DeviceType>();
  d_d2a = k_d2a.view<DeviceType>();
  d_d3a = k_d3a.view<DeviceType>();
  d_d4a = k_d4a.view<DeviceType>();
  d_zze = k_zze.view<DeviceType>();
  d_sw1 = k_sw1.view<DeviceType>();
  d_sw2 = k_sw2.view<DeviceType>();
  d_sw3 = k_sw3.view<DeviceType>();
  d_sw4 = k_sw4.view<DeviceType>();
  d_sw5 = k_sw5.view<DeviceType>();

  d_cutsq = typename AT::t_ffloat_2d_dl("pair_zbl:cutsq",n+1,n+1);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairZBLKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairZBL::init_one(i,j);

  k_z.h_view(i) = z[i];
  k_z.h_view(j) = z[j];
  k_d1a.h_view(i,j) = k_d1a.h_view(j,i) = d1a[i][j];
  k_d2a.h_view(i,j) = k_d2a.h_view(j,i) = d2a[i][j];
  k_d3a.h_view(i,j) = k_d3a.h_view(j,i) = d3a[i][j];
  k_d4a.h_view(i,j) = k_d4a.h_view(j,i) = d4a[i][j];
  k_zze.h_view(i,j) = k_zze.h_view(j,i) = zze[i][j];
  k_sw1.h_view(i,j) = k_sw1.h_view(j,i) = sw1[i][j];
  k_sw2.h_view(i,j) = k_sw2.h_view(j,i) = sw2[i][j];
  k_sw3.h_view(i,j) = k_sw3.h_view(j,i) = sw3[i][j];
  k_sw4.h_view(i,j) = k_sw4.h_view(j,i) = sw4[i][j];
  k_sw5.h_view(i,j) = k_sw5.h_view(j,i) = sw5[i][j];

  k_z.modify<LMPHostType>();
  k_d1a.modify<LMPHostType>();
  k_d2a.modify<LMPHostType>();
  k_d3a.modify<LMPHostType>();
  k_d4a.modify<LMPHostType>();
  k_zze.modify<LMPHostType>();
  k_sw1.modify<LMPHostType>();
  k_sw2.modify<LMPHostType>();
  k_sw3.modify<LMPHostType>();
  k_sw4.modify<LMPHostType>();
  k_sw5.modify<LMPHostType>();

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[i][j] = m_cutsq[j][i] = cutone*cutone;
  }

  return cutone;
}

/* ----------------------------------------------------------------------
   compute ZBL pair energy
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairZBLKokkos<DeviceType>::e_zbl(F_FLOAT r, int i, int j) const {

  const F_FLOAT d1aij = d_d1a(i,j);
  const F_FLOAT d2aij = d_d2a(i,j);
  const F_FLOAT d3aij = d_d3a(i,j);
  const F_FLOAT d4aij = d_d4a(i,j);
  const F_FLOAT zzeij = d_zze(i,j);
  const F_FLOAT rinv = 1.0/r;

  F_FLOAT sum = c1*exp(-d1aij*r);
  sum += c2*exp(-d2aij*r);
  sum += c3*exp(-d3aij*r);
  sum += c4*exp(-d4aij*r);

  F_FLOAT result = zzeij*sum*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL first derivative
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairZBLKokkos<DeviceType>::dzbldr(F_FLOAT r, int i, int j) const {

  const F_FLOAT d1aij = d_d1a(i,j);
  const F_FLOAT d2aij = d_d2a(i,j);
  const F_FLOAT d3aij = d_d3a(i,j);
  const F_FLOAT d4aij = d_d4a(i,j);
  const F_FLOAT zzeij = d_zze(i,j);
  const F_FLOAT rinv = 1.0/r;

  const F_FLOAT e1 = exp(-d1aij*r);
  const F_FLOAT e2 = exp(-d2aij*r);
  const F_FLOAT e3 = exp(-d3aij*r);
  const F_FLOAT e4 = exp(-d4aij*r);

  F_FLOAT sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  F_FLOAT sum_p = -c1*d1aij*e1;
  sum_p -= c2*d2aij*e2;
  sum_p -= c3*d3aij*e3;
  sum_p -= c4*d4aij*e4;

  F_FLOAT result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL second derivative
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairZBLKokkos<DeviceType>::d2zbldr2(F_FLOAT r, int i, int j) const {

  const F_FLOAT d1aij = d_d1a(i,j);
  const F_FLOAT d2aij = d_d2a(i,j);
  const F_FLOAT d3aij = d_d3a(i,j);
  const F_FLOAT d4aij = d_d4a(i,j);
  const F_FLOAT zzeij = d_zze(i,j);
  const F_FLOAT rinv = 1.0/r;

  const F_FLOAT e1 = exp(-d1aij*r);
  const F_FLOAT e2 = exp(-d2aij*r);
  const F_FLOAT e3 = exp(-d3aij*r);
  const F_FLOAT e4 = exp(-d4aij*r);

  F_FLOAT sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  F_FLOAT sum_p = c1*e1*d1aij;
  sum_p += c2*e2*d2aij;
  sum_p += c3*e3*d3aij;
  sum_p += c4*e4*d4aij;

  F_FLOAT sum_pp = c1*e1*d1aij*d1aij;
  sum_pp += c2*e2*d2aij*d2aij;
  sum_pp += c3*e3*d3aij*d3aij;
  sum_pp += c4*e4*d4aij*d4aij;

  F_FLOAT result = zzeij*(sum_pp + 2.0*sum_p*rinv +
                         2.0*sum*rinv*rinv)*rinv;

  return result;
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairZBLKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

namespace LAMMPS_NS {
template class PairZBLKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairZBLKokkos<LMPHostType>;
#endif
}
