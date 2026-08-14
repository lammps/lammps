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

#include "bond_mm3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondMM3Kokkos<DeviceType>::BondMM3Kokkos(LAMMPS *lmp) : BondMM3(lmp)
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
BondMM3Kokkos<DeviceType>::~BondMM3Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondMM3Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
    } else Kokkos::deep_copy(d_eatom,static_cast<KK_ACC_FLOAT>(0.0));
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"bond:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,static_cast<KK_ACC_FLOAT>(0.0));
  }

  k_r0.template sync<DeviceType>();
  k_k2.template sync<DeviceType>();

  // precompute K3 and K4 from host-side force->angstrom

  d_K3 = static_cast<KK_FLOAT>(-2.55 / force->angstrom);
  d_K4 = static_cast<KK_FLOAT>(7.0/12.0 * 2.55 * 2.55 /
    (force->angstrom * force->angstrom));

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagBondMM3Compute<1,1> >(0,nbondlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagBondMM3Compute<0,1> >(0,nbondlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagBondMM3Compute<1,0> >(0,nbondlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagBondMM3Compute<0,0> >(0,nbondlist),*this);
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
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondMM3Kokkos<DeviceType>::operator()(TagBondMM3Compute<NEWTON_BOND,EVFLAG>,
                                           const int &n, EV_FLOAT& ev) const
{
  const int i1 = bondlist(n,0);
  const int i2 = bondlist(n,1);
  const int type = bondlist(n,2);

  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  const KK_FLOAT r = sqrt(rsq);
  const KK_FLOAT dr = r - d_r0[type];
  const KK_FLOAT dr2 = dr*dr;

  const KK_FLOAT de_bond = static_cast<KK_FLOAT>(2.0)*d_k2[type]*dr*(
    static_cast<KK_FLOAT>(1.0)
    + static_cast<KK_FLOAT>(1.5)*d_K3*dr
    + static_cast<KK_FLOAT>(2.0)*d_K4*dr2);

  KK_FLOAT fbond = static_cast<KK_FLOAT>(0.0);
  if (r > static_cast<KK_FLOAT>(0.0)) fbond = -de_bond/r;

  KK_FLOAT ebond = static_cast<KK_FLOAT>(0.0);
  if (EVFLAG && eflag)
    ebond = d_k2[type]*dr2*(static_cast<KK_FLOAT>(1.0) + d_K3*dr + d_K4*dr2);

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
void BondMM3Kokkos<DeviceType>::operator()(TagBondMM3Compute<NEWTON_BOND,EVFLAG>,
                                           const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagBondMM3Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondMM3Kokkos<DeviceType>::allocate()
{
  BondMM3::allocate();
  int n = atom->nbondtypes;
  k_r0 = DAT::tdual_kkfloat_1d("BondMM3::r0",n+1);
  k_k2 = DAT::tdual_kkfloat_1d("BondMM3::k2",n+1);
  d_r0 = k_r0.template view<DeviceType>();
  d_k2 = k_k2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondMM3Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondMM3::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_r0.view_host()[i] = static_cast<KK_FLOAT>(r0[i]);
    k_k2.view_host()[i] = static_cast<KK_FLOAT>(k2[i]);
  }
  k_r0.modify_host();
  k_k2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondMM3Kokkos<DeviceType>::read_restart(FILE *fp)
{
  BondMM3::read_restart(fp);
  int n = atom->nbondtypes;
  for (int i = 1; i <= n; i++) {
    k_r0.view_host()[i] = static_cast<KK_FLOAT>(r0[i]);
    k_k2.view_host()[i] = static_cast<KK_FLOAT>(k2[i]);
  }
  k_r0.modify_host();
  k_k2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondMM3Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
    const KK_FLOAT &ebond, const KK_FLOAT &fbond,
    const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebond);
      else {
        const KK_ACC_FLOAT ebondhalf =
          static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
        if (i < nlocal) ev.evdwl += ebondhalf;
        if (j < nlocal) ev.evdwl += ebondhalf;
      }
    }
    if (eflag_atom) {
      const KK_ACC_FLOAT ebondhalf =
        static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
      if (newton_bond || i < nlocal) d_eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) d_eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v[6];
    v[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delx*fbond);
    v[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*dely*fbond);
    v[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delz*delz*fbond);
    v[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*dely*fbond);
    v[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delz*fbond);
    v[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*delz*fbond);

    if (vflag_global) {
      if (newton_bond) {
        for (int m = 0; m < 6; m++)
          ev.v[m] += static_cast<KK_ACC_FLOAT>(2.0)*v[m];
      } else {
        if (i < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
        if (j < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
      }
    }
    if (vflag_atom) {
      if (newton_bond || i < nlocal)
        for (int m = 0; m < 6; m++) d_vatom(i,m) += v[m];
      if (newton_bond || j < nlocal)
        for (int m = 0; m < 6; m++) d_vatom(j,m) += v[m];
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class BondMM3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondMM3Kokkos<LMPHostType>;
#endif
}
