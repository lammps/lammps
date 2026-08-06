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

#include "bond_gaussian_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>
#include <limits>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondGaussianKokkos<DeviceType>::BondGaussianKokkos(LAMMPS *lmp) : BondGaussian(lmp),
  kk_max_nterms(0)
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
BondGaussianKokkos<DeviceType>::~BondGaussianKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondGaussianKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  k_bond_temperature.template sync<DeviceType>();
  k_nterms.template sync<DeviceType>();
  k_alpha.template sync<DeviceType>();
  k_width.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondGaussianCompute<1,1> >(0,nbondlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondGaussianCompute<0,1> >(0,nbondlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondGaussianCompute<1,0> >(0,nbondlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondGaussianCompute<0,0> >(0,nbondlist),*this);
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
void BondGaussianKokkos<DeviceType>::operator()(TagBondGaussianCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // smallest normalized value of the accumulation type so the underflow guard
  // below works at any precision (casting the 2.3e-308 of the base class to
  // float gives 0.0)
  static constexpr KK_ACC_FLOAT SMALL_KK = std::numeric_limits<KK_ACC_FLOAT>::min();

  const int i1 = bondlist(n,0);
  const int i2 = bondlist(n,1);
  const int type = bondlist(n,2);

  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  const KK_FLOAT r = sqrt(rsq);

  const int nt = d_nterms[type];
  const KK_FLOAT kbT = d_bond_temperature[type];

  // the sums of the Gaussian terms must use the accumulation precision:
  // far out in the tails the individual terms underflow to zero in single
  // precision and log(sum_g_i) below would give -inf

  KK_ACC_FLOAT sum_g_i = static_cast<KK_ACC_FLOAT>(0.0);
  KK_ACC_FLOAT sum_numerator = static_cast<KK_ACC_FLOAT>(0.0);
  for (int i = 0; i < nt; i++) {
    const KK_ACC_FLOAT dr = r - d_r0(type,i);
    const KK_ACC_FLOAT wsq = static_cast<KK_ACC_FLOAT>(d_width(type,i)) * d_width(type,i);
    const KK_ACC_FLOAT prefactor = d_alpha(type,i) / (d_width(type,i) * sqrt(static_cast<KK_ACC_FLOAT>(MY_PI2)));
    const KK_ACC_FLOAT exponent = -static_cast<KK_ACC_FLOAT>(2.0) * dr * dr / wsq;
    const KK_ACC_FLOAT g_i = prefactor * exp(exponent);
    sum_g_i += g_i;
    sum_numerator += g_i * dr / wsq;
  }

  // avoid overflow
  if (sum_g_i < sum_numerator * SMALL_KK) sum_g_i = sum_numerator * SMALL_KK;

  // force & energy

  KK_FLOAT fbond = static_cast<KK_FLOAT>(0.0);
  if (r > static_cast<KK_FLOAT>(0.0))
    fbond = static_cast<KK_FLOAT>(-static_cast<KK_ACC_FLOAT>(4.0) * kbT * (sum_numerator / sum_g_i) / r);

  KK_FLOAT ebond = static_cast<KK_FLOAT>(0.0);
  if (eflag) ebond = static_cast<KK_FLOAT>(-kbT * log(sum_g_i));

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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondGaussianKokkos<DeviceType>::operator()(TagBondGaussianCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagBondGaussianCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondGaussianKokkos<DeviceType>::allocate()
{
  BondGaussian::allocate();

  int n = atom->nbondtypes;
  k_bond_temperature = DAT::tdual_kkfloat_1d("BondGaussian::bond_temperature",n+1);
  k_nterms = DAT::tdual_int_1d("BondGaussian::nterms",n+1);

  d_bond_temperature = k_bond_temperature.template view<DeviceType>();
  d_nterms = k_nterms.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondGaussianKokkos<DeviceType>::update_2d_views(int ntypes)
{
  // determine current max nterms across all bond types
  int new_max = 0;
  for (int i = 1; i <= ntypes; i++) {
    if (setflag[i] && nterms[i] > new_max) new_max = nterms[i];
  }
  if (new_max == 0) return;

  if (new_max > kk_max_nterms) {
    kk_max_nterms = new_max;
    k_alpha = DAT::tdual_kkfloat_2d("BondGaussian::alpha",ntypes+1,kk_max_nterms);
    k_width = DAT::tdual_kkfloat_2d("BondGaussian::width",ntypes+1,kk_max_nterms);
    k_r0 = DAT::tdual_kkfloat_2d("BondGaussian::r0",ntypes+1,kk_max_nterms);
    d_alpha = k_alpha.template view<DeviceType>();
    d_width = k_width.template view<DeviceType>();
    d_r0 = k_r0.template view<DeviceType>();
  }

  // populate host views with current data
  for (int i = 1; i <= ntypes; i++) {
    if (!setflag[i]) continue;
    // store kbT = k_B * T directly to avoid device-side access of force->boltz
    k_bond_temperature.view_host()[i] = static_cast<KK_FLOAT>(force->boltz * bond_temperature[i]);
    k_nterms.view_host()[i] = nterms[i];
    for (int j = 0; j < nterms[i]; j++) {
      k_alpha.view_host()(i,j) = static_cast<KK_FLOAT>(alpha[i][j]);
      k_width.view_host()(i,j) = static_cast<KK_FLOAT>(width[i][j]);
      k_r0.view_host()(i,j) = static_cast<KK_FLOAT>(r0[i][j]);
    }
  }

  k_bond_temperature.modify_host();
  k_nterms.modify_host();
  k_alpha.modify_host();
  k_width.modify_host();
  k_r0.modify_host();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void BondGaussianKokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondGaussian::coeff(narg, arg);

  update_2d_views(atom->nbondtypes);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void BondGaussianKokkos<DeviceType>::read_restart(FILE *fp)
{
  BondGaussian::read_restart(fp);

  update_2d_views(atom->nbondtypes);
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void BondGaussianKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
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
template class BondGaussianKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondGaussianKokkos<LMPHostType>;
#endif
}
