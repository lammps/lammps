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

#include "bond_table_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr int LINEAR_STYLE = 1;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondTableKokkos<DeviceType>::BondTableKokkos(LAMMPS *lmp) : BondTable(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_error_flag = DAT::tdual_int_scalar("BondTable::error_flag");
  d_error_flag = k_error_flag.template view<DeviceType>();
  h_error_flag = k_error_flag.view_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondTableKokkos<DeviceType>::~BondTableKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ----------------------------------------------------------------------
   copy the tabulated data built on the host by compute_table() into
   device views.  this runs once per setup, after every coeff() call, so
   restarts and repeated pair_coeff commands are covered
------------------------------------------------------------------------- */

template<class DeviceType>
void BondTableKokkos<DeviceType>::setup_tables()
{
  const int n = atom->nbondtypes;

  k_tabindex = DAT::tdual_int_1d("BondTable::tabindex",n+1);
  k_lo = DAT::tdual_kkfloat_1d("BondTable::lo",ntables);
  k_invdelta = DAT::tdual_kkfloat_1d("BondTable::invdelta",ntables);
  k_deltasq6 = DAT::tdual_kkfloat_1d("BondTable::deltasq6",ntables);

  // the arrays carry one extra element so that the spline branch can read
  // itable+1 at the last bin, where its weight is exactly zero

  k_r = DAT::tdual_kkfloat_2d("BondTable::r",ntables,tablength+1);
  k_e = DAT::tdual_kkfloat_2d("BondTable::e",ntables,tablength+1);
  k_de = DAT::tdual_kkfloat_2d("BondTable::de",ntables,tablength+1);
  k_f = DAT::tdual_kkfloat_2d("BondTable::f",ntables,tablength+1);
  k_df = DAT::tdual_kkfloat_2d("BondTable::df",ntables,tablength+1);
  k_e2 = DAT::tdual_kkfloat_2d("BondTable::e2",ntables,tablength+1);
  k_f2 = DAT::tdual_kkfloat_2d("BondTable::f2",ntables,tablength+1);

  for (int i = 1; i <= n; i++) k_tabindex.view_host()(i) = tabindex[i];

  for (int m = 0; m < ntables; m++) {
    const Table *tb = &tables[m];
    k_lo.view_host()(m) = static_cast<KK_FLOAT>(tb->lo);
    k_invdelta.view_host()(m) = static_cast<KK_FLOAT>(tb->invdelta);
    k_deltasq6.view_host()(m) = static_cast<KK_FLOAT>(tb->deltasq6);
    for (int i = 0; i < tablength; i++) {
      k_r.view_host()(m,i) = static_cast<KK_FLOAT>(tb->r[i]);
      k_e.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e[i]);
      k_de.view_host()(m,i) = static_cast<KK_FLOAT>(tb->de[i]);
      k_f.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f[i]);
      k_df.view_host()(m,i) = static_cast<KK_FLOAT>(tb->df[i]);
      k_e2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e2[i]);
      k_f2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f2[i]);
    }
  }

  k_tabindex.modify_host(); k_tabindex.template sync<DeviceType>();
  k_lo.modify_host(); k_lo.template sync<DeviceType>();
  k_invdelta.modify_host(); k_invdelta.template sync<DeviceType>();
  k_deltasq6.modify_host(); k_deltasq6.template sync<DeviceType>();
  k_r.modify_host(); k_r.template sync<DeviceType>();
  k_e.modify_host(); k_e.template sync<DeviceType>();
  k_de.modify_host(); k_de.template sync<DeviceType>();
  k_f.modify_host(); k_f.template sync<DeviceType>();
  k_df.modify_host(); k_df.template sync<DeviceType>();
  k_e2.modify_host(); k_e2.template sync<DeviceType>();
  k_f2.modify_host(); k_f2.template sync<DeviceType>();

  d_tabindex = k_tabindex.template view<DeviceType>();
  d_lo = k_lo.template view<DeviceType>();
  d_invdelta = k_invdelta.template view<DeviceType>();
  d_deltasq6 = k_deltasq6.template view<DeviceType>();
  d_r = k_r.template view<DeviceType>();
  d_e = k_e.template view<DeviceType>();
  d_de = k_de.template view<DeviceType>();
  d_f = k_f.template view<DeviceType>();
  d_df = k_df.template view<DeviceType>();
  d_e2 = k_e2.template view<DeviceType>();
  d_f2 = k_f2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondTableKokkos<DeviceType>::init_style()
{
  BondTable::init_style();

  setup_tables();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    if ((int)k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"bond:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_eatom,0.0);
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"bond:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,0.0);
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  h_error_flag() = 0;
  k_error_flag.modify_host();
  k_error_flag.template sync<DeviceType>();

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondTableCompute<1,1> >(0,nbondlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondTableCompute<0,1> >(0,nbondlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondTableCompute<1,0> >(0,nbondlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondTableCompute<0,0> >(0,nbondlist),*this);
    }
  }

  // error check

  k_error_flag.template modify<DeviceType>();
  k_error_flag.sync_host();
  if (h_error_flag() == 1)
    error->one(FLERR,"Bond length < table inner cutoff");
  else if (h_error_flag() == 2)
    error->one(FLERR,"Bond length > table outer cutoff");

  if (eflag_global) energy += static_cast<double>(ev.evdwl);
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

  copymode = 0;
}

/* ----------------------------------------------------------------------
   device version of BondTable::uf_lookup()
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondTableKokkos<DeviceType>::uf_lookup_kk(const int &type, const KK_FLOAT &x_in,
                                               KK_FLOAT &u, KK_FLOAT &mdu) const
{
  u = mdu = static_cast<KK_FLOAT>(0.0);

  const int tb = d_tabindex[type];
  const KK_FLOAT lo = d_lo[tb];
  const KK_FLOAT invdelta = d_invdelta[tb];
  const int itable = static_cast<int>((x_in - lo) * invdelta);

  if (itable < 0) { d_error_flag() = 1; return; }
  if (itable >= tablength) { d_error_flag() = 2; return; }

  const KK_FLOAT b = (x_in - d_r(tb,itable)) * invdelta;

  if (tabstyle == LINEAR_STYLE) {
    u = d_e(tb,itable) + b * d_de(tb,itable);
    mdu = d_f(tb,itable) + b * d_df(tb,itable);
  } else {
    const KK_FLOAT a = static_cast<KK_FLOAT>(1.0) - b;
    const KK_FLOAT deltasq6 = d_deltasq6[tb];
    u = a * d_e(tb,itable) + b * d_e(tb,itable+1) +
        ((a*a*a - a) * d_e2(tb,itable) + (b*b*b - b) * d_e2(tb,itable+1)) * deltasq6;
    mdu = a * d_f(tb,itable) + b * d_f(tb,itable+1) +
        ((a*a*a - a) * d_f2(tb,itable) + (b*b*b - b) * d_f2(tb,itable+1)) * deltasq6;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondTableKokkos<DeviceType>::operator()(TagBondTableCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  const int i1 = bondlist(n,0);
  const int i2 = bondlist(n,1);
  const int type = bondlist(n,2);

  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  const KK_FLOAT r = Kokkos::sqrt(rsq);

  // force & energy

  KK_FLOAT u,mdu;
  uf_lookup_kk(type,r,u,mdu);

  const KK_FLOAT fbond = mdu / r;
  const KK_FLOAT ebond = u;

  // apply force to each of 2 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += static_cast<KK_ACC_FLOAT>(delx*fbond);
    f(i1,1) += static_cast<KK_ACC_FLOAT>(dely*fbond);
    f(i1,2) += static_cast<KK_ACC_FLOAT>(delz*fbond);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) -= static_cast<KK_ACC_FLOAT>(delx*fbond);
    f(i2,1) -= static_cast<KK_ACC_FLOAT>(dely*fbond);
    f(i2,2) -= static_cast<KK_ACC_FLOAT>(delz*fbond);
  }

  if (EVFLAG) ev_tally(ev,i1,i2,ebond,fbond,delx,dely,delz);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondTableKokkos<DeviceType>::operator()(TagBondTableCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagBondTableCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondTableKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &ebond, const KK_FLOAT &fbond, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebond);
      else {
        KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
        if (i < nlocal) ev.evdwl += ebondhalf;
        if (j < nlocal) ev.evdwl += ebondhalf;
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
      if (newton_bond || i < nlocal) d_eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) d_eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v_half_acc[6];
    v_half_acc[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delx*fbond);
    v_half_acc[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*dely*fbond);
    v_half_acc[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delz*delz*fbond);
    v_half_acc[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*dely*fbond);
    v_half_acc[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delz*fbond);
    v_half_acc[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*delz*fbond);

    if (vflag_global) {
      if (newton_bond) {
        for (int k = 0; k < 6; k++) ev.v[k] += static_cast<KK_ACC_FLOAT>(2.0)*v_half_acc[k];
      } else {
        if (i < nlocal) for (int k = 0; k < 6; k++) ev.v[k] += v_half_acc[k];
        if (j < nlocal) for (int k = 0; k < 6; k++) ev.v[k] += v_half_acc[k];
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal)
        for (int k = 0; k < 6; k++) d_vatom(i,k) += v_half_acc[k];
      if (newton_bond || j < nlocal)
        for (int k = 0; k < 6; k++) d_vatom(j,k) += v_half_acc[k];
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class BondTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondTableKokkos<LMPHostType>;
#endif
}
