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

#include "pair_zbl_kokkos.h"
#include <cmath>
#include <cstring>
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
#include "kokkos.h"

// From J.F. Zeigler, J. P. Biersack and U. Littmark,
// "The Stopping and Range of Ions in Matter" volume 1, Pergamon, 1985.

using namespace LAMMPS_NS;
using namespace PairZBLConstants;

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairZBLKokkos<Space>::PairZBLKokkos(LAMMPS *lmp) : PairZBL(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairZBLKokkos<Space>::~PairZBLKokkos()
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

template<ExecutionSpace Space>
void PairZBLKokkos<Space>::init_style()
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
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }

  Kokkos::deep_copy(d_cutsq,cut_globalsq);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairZBLKokkos<Space>::compute(int eflag_in, int vflag_in)
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
  if (eflag || vflag) atomKK->modified(Space,datamask_modify);
  else atomKK->modified(Space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  DualViewHelper<Space>::sync(k_z);
  DualViewHelper<Space>::sync(k_d1a);
  DualViewHelper<Space>::sync(k_d2a);
  DualViewHelper<Space>::sync(k_d3a);
  DualViewHelper<Space>::sync(k_d4a);
  DualViewHelper<Space>::sync(k_zze);
  DualViewHelper<Space>::sync(k_sw1);
  DualViewHelper<Space>::sync(k_sw2);
  DualViewHelper<Space>::sync(k_sw3);
  DualViewHelper<Space>::sync(k_sw4);
  DualViewHelper<Space>::sync(k_sw5);

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<Space,PairZBLKokkos<Space>,void >(this,(NeighListKokkos<Space>*)list);

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
KK_FLOAT PairZBLKokkos<Space>::
compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const KK_FLOAT r = sqrt(rsq);
  KK_FLOAT fpair = dzbldr(r, itype, jtype);

  if (rsq > cut_innersq) {
    const KK_FLOAT t = r - cut_inner;
    const KK_FLOAT fswitch = t*t *
           (d_sw1(itype,jtype) + d_sw2(itype,jtype)*t);
    fpair += fswitch;
  }

  fpair *= -1.0/r;
  return fpair;
}

template<ExecutionSpace Space>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<Space>::
compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const KK_FLOAT r = sqrt(rsq);
  KK_FLOAT evdwl = e_zbl(r, itype, jtype);
  evdwl += d_sw5(itype,jtype);
  if (rsq > cut_innersq) {
    const KK_FLOAT t = r - cut_inner;
    const KK_FLOAT eswitch = t*t*t *
      (d_sw3(itype,jtype) + d_sw4(itype,jtype)*t);
    evdwl += eswitch;
  }
  return evdwl;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairZBLKokkos<Space>::allocate()
{
  PairZBL::allocate();

  int n = atom->ntypes;

  k_z   = DAT::tdual_float_1d("pair_zbl:z  ",n+1);
  k_d1a = DAT::tdual_float_2d_dl("pair_zbl:d1a",n+1,n+1);
  k_d2a = DAT::tdual_float_2d_dl("pair_zbl:d2a",n+1,n+1);
  k_d3a = DAT::tdual_float_2d_dl("pair_zbl:d3a",n+1,n+1);
  k_d4a = DAT::tdual_float_2d_dl("pair_zbl:d4a",n+1,n+1);
  k_zze = DAT::tdual_float_2d_dl("pair_zbl:zze",n+1,n+1);
  k_sw1 = DAT::tdual_float_2d_dl("pair_zbl:sw1",n+1,n+1);
  k_sw2 = DAT::tdual_float_2d_dl("pair_zbl:sw2",n+1,n+1);
  k_sw3 = DAT::tdual_float_2d_dl("pair_zbl:sw3",n+1,n+1);
  k_sw4 = DAT::tdual_float_2d_dl("pair_zbl:sw4",n+1,n+1);
  k_sw5 = DAT::tdual_float_2d_dl("pair_zbl:sw5",n+1,n+1);

  d_z   = DualViewHelper<Space>::view(k_z);
  d_d1a = DualViewHelper<Space>::view(k_d1a);
  d_d2a = DualViewHelper<Space>::view(k_d2a);
  d_d3a = DualViewHelper<Space>::view(k_d3a);
  d_d4a = DualViewHelper<Space>::view(k_d4a);
  d_zze = DualViewHelper<Space>::view(k_zze);
  d_sw1 = DualViewHelper<Space>::view(k_sw1);
  d_sw2 = DualViewHelper<Space>::view(k_sw2);
  d_sw3 = DualViewHelper<Space>::view(k_sw3);
  d_sw4 = DualViewHelper<Space>::view(k_sw4);
  d_sw5 = DualViewHelper<Space>::view(k_sw5);

  d_cutsq = typename AT::t_float_2d_dl("pair_zbl:cutsq",n+1,n+1);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairZBLKokkos<Space>::init_one(int i, int j)
{
  KK_FLOAT cutone = PairZBL::init_one(i,j);

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

  k_z.modify_host();
  k_d1a.modify_host();
  k_d2a.modify_host();
  k_d3a.modify_host();
  k_d4a.modify_host();
  k_zze.modify_host();
  k_sw1.modify_host();
  k_sw2.modify_host();
  k_sw3.modify_host();
  k_sw4.modify_host();
  k_sw5.modify_host();

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[i][j] = m_cutsq[j][i] = cutone*cutone;
  }

  return cutone;
}

/* ----------------------------------------------------------------------
   compute ZBL pair energy
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<Space>::e_zbl(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = 1.0/r;

  KK_FLOAT sum = c1*exp(-d1aij*r);
  sum += c2*exp(-d2aij*r);
  sum += c3*exp(-d3aij*r);
  sum += c4*exp(-d4aij*r);

  KK_FLOAT result = zzeij*sum*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL first derivative
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<Space>::dzbldr(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = 1.0/r;

  const KK_FLOAT e1 = exp(-d1aij*r);
  const KK_FLOAT e2 = exp(-d2aij*r);
  const KK_FLOAT e3 = exp(-d3aij*r);
  const KK_FLOAT e4 = exp(-d4aij*r);

  KK_FLOAT sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  KK_FLOAT sum_p = -c1*d1aij*e1;
  sum_p -= c2*d2aij*e2;
  sum_p -= c3*d3aij*e3;
  sum_p -= c4*d4aij*e4;

  KK_FLOAT result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL second derivative
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<Space>::d2zbldr2(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = 1.0/r;

  const KK_FLOAT e1 = exp(-d1aij*r);
  const KK_FLOAT e2 = exp(-d2aij*r);
  const KK_FLOAT e3 = exp(-d3aij*r);
  const KK_FLOAT e4 = exp(-d4aij*r);

  KK_FLOAT sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  KK_FLOAT sum_p = c1*e1*d1aij;
  sum_p += c2*e2*d2aij;
  sum_p += c3*e3*d3aij;
  sum_p += c4*e4*d4aij;

  KK_FLOAT sum_pp = c1*e1*d1aij*d1aij;
  sum_pp += c2*e2*d2aij*d2aij;
  sum_pp += c3*e3*d3aij*d3aij;
  sum_pp += c4*e4*d4aij*d4aij;

  KK_FLOAT result = zzeij*(sum_pp + 2.0*sum_p*rinv +
                         2.0*sum*rinv*rinv)*rinv;

  return result;
}


/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairZBLKokkos<Space>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

namespace LAMMPS_NS {
template class PairZBLKokkos<Device>;
template class PairZBLKokkos<Host>;
}
