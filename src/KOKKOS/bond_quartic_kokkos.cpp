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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "bond_quartic_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_CUBEROOT2;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondQuarticKokkos<DeviceType>::BondQuarticKokkos(LAMMPS *lmp) : BondQuartic(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondQuarticKokkos<DeviceType>::~BondQuarticKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondQuarticKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  k_k.template sync<DeviceType>();
  k_b1.template sync<DeviceType>();
  k_b2.template sync<DeviceType>();
  k_rc.template sync<DeviceType>();
  k_u0.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  // allocate/resize brokenflag array

  if ((int)k_brokenflag.extent(0) < nbondlist) {
    k_brokenflag = DAT::tdual_int_1d("BondQuartic::brokenflag",nbondlist);
    d_brokenflag = k_brokenflag.template view<DeviceType>();
  }
  Kokkos::deep_copy(d_brokenflag,0);

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondQuarticCompute<1,1> >(0,nbondlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondQuarticCompute<0,1> >(0,nbondlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondQuarticCompute<1,0> >(0,nbondlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondQuarticCompute<0,0> >(0,nbondlist),*this);
    }
  }

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

  // mark forces as modified by the kernels above, then pull positions,
  // forces, and bond topology into the legacy host arrays for the bond
  // breaking and pair correction post-processing below

  atomKK->modified(execution_space, F_MASK);
  atomKK->sync(Host, X_MASK | F_MASK | BOND_MASK | TAG_MASK | TYPE_MASK);
  atomKK->k_f.modify_host();

  k_brokenflag.sync_host();

  // host-side post-processing: bond breaking and pair correction

  neighborKK->k_bondlist.sync_host();
  auto h_bondlist = neighborKK->k_bondlist.view_host();
  auto h_brokenflag = k_brokenflag.view_host();

  double **x_host = atom->x;
  double **f_host = atom->f;

  // ensure pair->ev_tally() will use 1-4 virial contribution

  if (vflag_global == VIRIAL_FDOTR)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  double **cutsq = force->pair->cutsq;

  for (int n = 0; n < nbondlist; n++) {

    // skip already-broken bonds (type <= 0) and unflagged bonds
    if (h_bondlist(n,2) <= 0) continue;

    if (h_brokenflag(n)) {
      int i1 = h_bondlist(n,0);
      int i2 = h_bondlist(n,1);
      h_bondlist(n,2) = 0;
      for (int m = 0; m < atom->num_bond[i1]; m++)
        if (atom->bond_atom[i1][m] == atom->tag[i2]) atom->bond_type[i1][m] = 0;
      if (i2 < atom->nlocal)
        for (int m = 0; m < atom->num_bond[i2]; m++)
          if (atom->bond_atom[i2][m] == atom->tag[i1]) atom->bond_type[i2][m] = 0;
      continue;
    }

    int i1 = h_bondlist(n,0);
    int i2 = h_bondlist(n,1);
    int type_n = h_bondlist(n,2);

    double delx = x_host[i1][0] - x_host[i2][0];
    double dely = x_host[i1][1] - x_host[i2][1];
    double delz = x_host[i1][2] - x_host[i2][2];
    double rsq = delx*delx + dely*dely + delz*delz;

    int itype = atom->type[i1];
    int jtype = atom->type[i2];

    if (rsq < cutsq[itype][jtype]) {
      double evdwl = 0.0, fpair = 0.0;
      evdwl = -force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,fpair);
      fpair = -fpair;

      if (newton_bond || i1 < atom->nlocal) {
        f_host[i1][0] += delx*fpair;
        f_host[i1][1] += dely*fpair;
        f_host[i1][2] += delz*fpair;
      }
      if (newton_bond || i2 < atom->nlocal) {
        f_host[i2][0] -= delx*fpair;
        f_host[i2][1] -= dely*fpair;
        f_host[i2][2] -= delz*fpair;
      }

      if (evflag)
        force->pair->ev_tally(i1,i2,atom->nlocal,newton_bond,evdwl,0.0,fpair,delx,dely,delz);
    }

    // suppress unused variable warning
    (void) type_n;
  }

  // mark bondlist and broken bond topology as modified on host

  neighborKK->k_bondlist.modify_host();
  atomKK->modified(Host, BOND_MASK);

  // transform the updated forces back into the Kokkos views, since the
  // integrator marks them as modified in execution_space after compute()

  atomKK->sync(execution_space, F_MASK);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondQuarticKokkos<DeviceType>::operator()(TagBondQuarticCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // skip already-broken bonds

  if (bondlist(n,2) <= 0) return;

  const int i1 = bondlist(n,0);
  const int i2 = bondlist(n,1);
  const int type = bondlist(n,2);

  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

  // if bond breaks, flag it and return (no force)

  if (rsq > d_rc[type]*d_rc[type]) {
    d_brokenflag(n) = 1;
    return;
  }

  const KK_FLOAT r = sqrt(rsq);
  const KK_FLOAT dr = r - d_rc[type];
  const KK_FLOAT r2 = dr*dr;
  const KK_FLOAT ra = dr - d_b1[type];
  const KK_FLOAT rb = dr - d_b2[type];

  KK_FLOAT fbond = -d_k[type]/r *
      (r2*(ra+rb) + static_cast<KK_FLOAT>(2.0)*dr*ra*rb);

  KK_FLOAT sr6 = static_cast<KK_FLOAT>(0.0);
  if (rsq < static_cast<KK_FLOAT>(MY_CUBEROOT2)) {
    const KK_FLOAT sr2 = static_cast<KK_FLOAT>(1.0)/rsq;
    sr6 = sr2*sr2*sr2;
    fbond += static_cast<KK_FLOAT>(48.0)*sr6*(sr6 - static_cast<KK_FLOAT>(0.5))/rsq;
  }

  KK_FLOAT ebond = static_cast<KK_FLOAT>(0.0);
  if (eflag) {
    ebond = d_k[type]*r2*ra*rb + d_u0[type];
    if (rsq < static_cast<KK_FLOAT>(MY_CUBEROOT2))
      ebond += static_cast<KK_FLOAT>(4.0)*sr6*(sr6 - static_cast<KK_FLOAT>(1.0)) +
          static_cast<KK_FLOAT>(1.0);
  }

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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondQuarticKokkos<DeviceType>::operator()(TagBondQuarticCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagBondQuarticCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondQuarticKokkos<DeviceType>::allocate()
{
  BondQuartic::allocate();

  int n = atom->nbondtypes;
  k_k = DAT::tdual_kkfloat_1d("BondQuartic::k",n+1);
  k_b1 = DAT::tdual_kkfloat_1d("BondQuartic::b1",n+1);
  k_b2 = DAT::tdual_kkfloat_1d("BondQuartic::b2",n+1);
  k_rc = DAT::tdual_kkfloat_1d("BondQuartic::rc",n+1);
  k_u0 = DAT::tdual_kkfloat_1d("BondQuartic::u0",n+1);

  d_k = k_k.template view<DeviceType>();
  d_b1 = k_b1.template view<DeviceType>();
  d_b2 = k_b2.template view<DeviceType>();
  d_rc = k_rc.template view<DeviceType>();
  d_u0 = k_u0.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void BondQuarticKokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondQuartic::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_b1.view_host()[i] = static_cast<KK_FLOAT>(b1[i]);
    k_b2.view_host()[i] = static_cast<KK_FLOAT>(b2[i]);
    k_rc.view_host()[i] = static_cast<KK_FLOAT>(rc[i]);
    k_u0.view_host()[i] = static_cast<KK_FLOAT>(u0[i]);
  }

  k_k.modify_host();
  k_b1.modify_host();
  k_b2.modify_host();
  k_rc.modify_host();
  k_u0.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void BondQuarticKokkos<DeviceType>::read_restart(FILE *fp)
{
  BondQuartic::read_restart(fp);

  int n = atom->nbondtypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_b1.view_host()[i] = static_cast<KK_FLOAT>(b1[i]);
    k_b2.view_host()[i] = static_cast<KK_FLOAT>(b2[i]);
    k_rc.view_host()[i] = static_cast<KK_FLOAT>(rc[i]);
    k_u0.view_host()[i] = static_cast<KK_FLOAT>(u0[i]);
  }

  k_k.modify_host();
  k_b1.modify_host();
  k_b2.modify_host();
  k_rc.modify_host();
  k_u0.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondQuarticKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &ebond, const KK_FLOAT &fbond, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebond);
      else {
        KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
        if (i < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebondhalf);
        if (j < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebondhalf);
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
      if (newton_bond || i < nlocal) d_eatom[i] += static_cast<KK_ACC_FLOAT>(ebondhalf);
      if (newton_bond || j < nlocal) d_eatom[j] += static_cast<KK_ACC_FLOAT>(ebondhalf);
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
        for (int n = 0; n < 6; n++)
          ev.v[n] += static_cast<KK_ACC_FLOAT>(2.0)*v_half_acc[n];
      } else {
        if (i < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_half_acc[n];
        }
        if (j < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_half_acc[n];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i,n) += v_half_acc[n];
      }
      if (newton_bond || j < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(j,n) += v_half_acc[n];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class BondQuarticKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondQuarticKokkos<LMPHostType>;
#endif
}
