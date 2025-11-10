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
   Contributing authors: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_zbl_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include "pair_zbl_const.h"

#include <cmath>
#include <cstring>

// From J.F. Zeigler, J. P. Biersack and U. Littmark,
// "The Stopping and Range of Ions in Matter" volume 1, Pergamon, 1985.

using namespace LAMMPS_NS;
using namespace PairZBLConstants;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairZBLKokkos<DeviceType>::PairZBLKokkos(LAMMPS *lmp) : PairZBL(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
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
  }
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairZBLKokkos<DeviceType>::init_style()
{
  PairZBL::init_style();

  Kokkos::deep_copy(d_cutsq,static_cast<KK_FLOAT>(cut_globalsq));

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
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
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  c1_kk = static_cast<KK_FLOAT>(c1);
  c2_kk = static_cast<KK_FLOAT>(c2);
  c3_kk = static_cast<KK_FLOAT>(c3);
  c4_kk = static_cast<KK_FLOAT>(c4);

  cut_inner_kk = static_cast<KK_FLOAT>(cut_inner);
  cut_innersq_kk = static_cast<KK_FLOAT>(cut_innersq);

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

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = sqrt(rsq);
  KK_FLOAT fpair = dzbldr(r, itype, jtype);

  if (rsq > cut_innersq_kk) {
    const KK_FLOAT t = r - cut_inner_kk;
    const KK_FLOAT fswitch = t*t *
           (d_sw1(itype,jtype) + d_sw2(itype,jtype)*t);
    fpair += fswitch;
  }

  fpair *= -static_cast<KK_FLOAT>(1.0) / r;
  return fpair;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = sqrt(rsq);
  KK_FLOAT evdwl = e_zbl(r, itype, jtype);
  evdwl += d_sw5(itype,jtype);
  if (rsq > cut_innersq_kk) {
    const KK_FLOAT t = r - cut_inner_kk;
    const KK_FLOAT eswitch = t*t*t *
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

  k_z   = DAT::tdual_kkfloat_1d("pair_zbl:z  ",n+1);
  k_d1a = DAT::tdual_kkfloat_2d_dl("pair_zbl:d1a",n+1,n+1);
  k_d2a = DAT::tdual_kkfloat_2d_dl("pair_zbl:d2a",n+1,n+1);
  k_d3a = DAT::tdual_kkfloat_2d_dl("pair_zbl:d3a",n+1,n+1);
  k_d4a = DAT::tdual_kkfloat_2d_dl("pair_zbl:d4a",n+1,n+1);
  k_zze = DAT::tdual_kkfloat_2d_dl("pair_zbl:zze",n+1,n+1);
  k_sw1 = DAT::tdual_kkfloat_2d_dl("pair_zbl:sw1",n+1,n+1);
  k_sw2 = DAT::tdual_kkfloat_2d_dl("pair_zbl:sw2",n+1,n+1);
  k_sw3 = DAT::tdual_kkfloat_2d_dl("pair_zbl:sw3",n+1,n+1);
  k_sw4 = DAT::tdual_kkfloat_2d_dl("pair_zbl:sw4",n+1,n+1);
  k_sw5 = DAT::tdual_kkfloat_2d_dl("pair_zbl:sw5",n+1,n+1);

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

  d_cutsq = typename AT::t_kkfloat_2d_dl("pair_zbl:cutsq",n+1,n+1);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairZBLKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairZBL::init_one(i,j);

  k_z.view_host()(i) = static_cast<KK_FLOAT>(z[i]);
  k_z.view_host()(j) = static_cast<KK_FLOAT>(z[j]);
  k_d1a.view_host()(i,j) = k_d1a.view_host()(j,i) = static_cast<KK_FLOAT>(d1a[i][j]);
  k_d2a.view_host()(i,j) = k_d2a.view_host()(j,i) = static_cast<KK_FLOAT>(d2a[i][j]);
  k_d3a.view_host()(i,j) = k_d3a.view_host()(j,i) = static_cast<KK_FLOAT>(d3a[i][j]);
  k_d4a.view_host()(i,j) = k_d4a.view_host()(j,i) = static_cast<KK_FLOAT>(d4a[i][j]);
  k_zze.view_host()(i,j) = k_zze.view_host()(j,i) = static_cast<KK_FLOAT>(zze[i][j]);
  k_sw1.view_host()(i,j) = k_sw1.view_host()(j,i) = static_cast<KK_FLOAT>(sw1[i][j]);
  k_sw2.view_host()(i,j) = k_sw2.view_host()(j,i) = static_cast<KK_FLOAT>(sw2[i][j]);
  k_sw3.view_host()(i,j) = k_sw3.view_host()(j,i) = static_cast<KK_FLOAT>(sw3[i][j]);
  k_sw4.view_host()(i,j) = k_sw4.view_host()(j,i) = static_cast<KK_FLOAT>(sw4[i][j]);
  k_sw5.view_host()(i,j) = k_sw5.view_host()(j,i) = static_cast<KK_FLOAT>(sw5[i][j]);

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

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[i][j] = m_cutsq[j][i] = static_cast<KK_FLOAT>(cutone*cutone);
  }

  return cutone;
}

/* ----------------------------------------------------------------------
   compute ZBL pair energy
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<DeviceType>::e_zbl(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;

  KK_FLOAT sum = c1_kk*exp(-d1aij*r);
  sum += c2_kk*exp(-d2aij*r);
  sum += c3_kk*exp(-d3aij*r);
  sum += c4_kk*exp(-d4aij*r);

  KK_FLOAT result = zzeij*sum*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL first derivative
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<DeviceType>::dzbldr(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;

  const KK_FLOAT e1 = exp(-d1aij*r);
  const KK_FLOAT e2 = exp(-d2aij*r);
  const KK_FLOAT e3 = exp(-d3aij*r);
  const KK_FLOAT e4 = exp(-d4aij*r);

  KK_FLOAT sum = c1_kk*e1;
  sum += c2_kk*e2;
  sum += c3_kk*e3;
  sum += c4_kk*e4;

  KK_FLOAT sum_p = -c1_kk*d1aij*e1;
  sum_p -= c2_kk*d2aij*e2;
  sum_p -= c3_kk*d3aij*e3;
  sum_p -= c4_kk*d4aij*e4;

  KK_FLOAT result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL second derivative
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairZBLKokkos<DeviceType>::d2zbldr2(KK_FLOAT r, int i, int j) const {

  const KK_FLOAT d1aij = d_d1a(i,j);
  const KK_FLOAT d2aij = d_d2a(i,j);
  const KK_FLOAT d3aij = d_d3a(i,j);
  const KK_FLOAT d4aij = d_d4a(i,j);
  const KK_FLOAT zzeij = d_zze(i,j);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;

  const KK_FLOAT e1 = exp(-d1aij*r);
  const KK_FLOAT e2 = exp(-d2aij*r);
  const KK_FLOAT e3 = exp(-d3aij*r);
  const KK_FLOAT e4 = exp(-d4aij*r);

  KK_FLOAT sum = c1_kk*e1;
  sum += c2_kk*e2;
  sum += c3_kk*e3;
  sum += c4_kk*e4;

  KK_FLOAT sum_p = c1_kk*e1*d1aij;
  sum_p += c2_kk*e2*d2aij;
  sum_p += c3_kk*e3*d3aij;
  sum_p += c4_kk*e4*d4aij;

  KK_FLOAT sum_pp = c1_kk*e1*d1aij*d1aij;
  sum_pp += c2_kk*e2*d2aij*d2aij;
  sum_pp += c3_kk*e3*d3aij*d3aij;
  sum_pp += c4_kk*e4*d4aij*d4aij;

  KK_FLOAT result = zzeij*(sum_pp + static_cast<KK_FLOAT>(2.0)*sum_p*rinv +
                         static_cast<KK_FLOAT>(2.0)*sum*rinv*rinv)*rinv;

  return result;
}

namespace LAMMPS_NS {
template class PairZBLKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairZBLKokkos<LMPHostType>;
#endif
}
