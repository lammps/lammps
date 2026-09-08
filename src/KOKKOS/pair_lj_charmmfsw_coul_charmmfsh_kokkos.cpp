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
   Contributing author: Mitch Murphy (alphataubio@gmail.com)

   Based on serial kspace lj-fsw sections (force-switched) provided by
   Robert Meissner and Lucio Colombi Ciacchi of Bremen University, Germany,
   with additional assistance from Robert A. Latour, Clemson University

 ------------------------------------------------------------------------- */

#include "pair_lj_charmmfsw_coul_charmmfsh_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::PairLJCharmmfswCoulCharmmfshKokkos(LAMMPS *lmp):PairLJCharmmfswCoulCharmmfsh(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::~PairLJCharmmfswCoulCharmmfshKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cutsq.template sync<DeviceType>();
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
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  ev = pair_compute<PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>,CoulLongTable<0> >
      (this,(NeighListKokkos<DeviceType>*)list);


  if (eflag) {
    eng_vdwl += static_cast<double>(ev.evdwl);
    eng_coul += static_cast<double>(ev.ecoul);
  }
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

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const KK_FLOAT cut_lj_innersq_kk = static_cast<KK_FLOAT>(cut_lj_innersq);
  const KK_FLOAT cut_ljsq_kk = static_cast<KK_FLOAT>(cut_ljsq);
  const KK_FLOAT denom_lj_kk = static_cast<KK_FLOAT>(denom_lj);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  KK_FLOAT forcelj, switch1;

  forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  if (rsq > cut_lj_innersq_kk) {
    switch1 = (cut_ljsq_kk-rsq) * (cut_ljsq_kk-rsq) *
              (cut_ljsq_kk + static_cast<KK_FLOAT>(2.0)*rsq - static_cast<KK_FLOAT>(3.0)*cut_lj_innersq_kk) / denom_lj_kk;
    forcelj = forcelj*switch1;
  }

  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const KK_FLOAT cut_lj_innersq_kk = static_cast<KK_FLOAT>(cut_lj_innersq);
  const KK_FLOAT cut_lj6_kk = static_cast<KK_FLOAT>(cut_lj6);
  const KK_FLOAT cut_lj6inv_kk = static_cast<KK_FLOAT>(cut_lj6inv);
  const KK_FLOAT cut_lj3_kk = static_cast<KK_FLOAT>(cut_lj3);
  const KK_FLOAT cut_lj3inv_kk = static_cast<KK_FLOAT>(cut_lj3inv);
  const KK_FLOAT cut_lj_inner6inv_kk = static_cast<KK_FLOAT>(cut_lj_inner6inv);
  const KK_FLOAT cut_lj_inner3inv_kk = static_cast<KK_FLOAT>(cut_lj_inner3inv);
  const KK_FLOAT denom_lj12_kk = static_cast<KK_FLOAT>(denom_lj12);
  const KK_FLOAT denom_lj6_kk = static_cast<KK_FLOAT>(denom_lj6);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0)/r;
  const KK_FLOAT r3inv = rinv*rinv*rinv;
  KK_FLOAT englj, englj12, englj6;

  if (rsq > cut_lj_innersq_kk) {
    englj12 = (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*cut_lj6_kk*
      denom_lj12_kk * (r6inv - cut_lj6inv_kk)*(r6inv - cut_lj6inv_kk);
    englj6 = -(STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4)*
      cut_lj3_kk*denom_lj6_kk * (r3inv - cut_lj3inv_kk)*(r3inv - cut_lj3inv_kk);
    englj = englj12 + englj6;
  } else {
    englj12 = r6inv*(STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
    (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*cut_lj_inner6inv_kk*cut_lj6inv_kk;
    englj6 = -(STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4)*r6inv +
      (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4)*
      cut_lj_inner3inv_kk*cut_lj3inv_kk;
    englj = englj12 + englj6;
  }
  return englj;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS,  class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT cut_coulinv_kk = static_cast<KK_FLOAT>(cut_coulinv);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r = Kokkos::sqrt(rsq);

  const KK_FLOAT forcecoul = qqrd2e * qtmp*q(j) *
    (Kokkos::sqrt(r2inv) - r*cut_coulinv_kk*cut_coulinv_kk);

  return factor_coul*forcecoul*r2inv;
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT cut_coulinv_kk = static_cast<KK_FLOAT>(cut_coulinv);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r = Kokkos::sqrt(rsq);

  return factor_coul * qqrd2e * qtmp*q(j) *
    (Kokkos::sqrt(r2inv) + cut_coulinv_kk*cut_coulinv_kk*r -
     static_cast<KK_FLOAT>(2.0)*cut_coulinv_kk);
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::allocate()
{
  PairLJCharmmfswCoulCharmmfsh::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq = typename AT::t_kkfloat_2d("pair:cut_ljsq",n+1,n+1);

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType>("PairLJCharmmfswCoulCharmmfsh::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::init_style()
{
  PairLJCharmmfswCoulCharmmfsh::init_style();

  Kokkos::deep_copy(d_cut_ljsq,static_cast<KK_FLOAT>(cut_ljsq));
  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

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

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCharmmfswCoulCharmmfshKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCharmmfswCoulCharmmfsh::init_one(i,j);

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsq);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsq);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCharmmfswCoulCharmmfshKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCharmmfswCoulCharmmfshKokkos<LMPHostType>;
#endif
}
