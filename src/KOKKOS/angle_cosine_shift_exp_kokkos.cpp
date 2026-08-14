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

#include "angle_cosine_shift_exp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleCosineShiftExpKokkos<DeviceType>::AngleCosineShiftExpKokkos(LAMMPS *lmp) : AngleCosineShiftExp(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleCosineShiftExpKokkos<DeviceType>::~AngleCosineShiftExpKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCosineShiftExpKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    if ((int)k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_eatom,0.0);
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"angle:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,0.0);
  }

  k_umin.template sync<DeviceType>();
  k_a.template sync<DeviceType>();
  k_opt1.template sync<DeviceType>();
  k_sint.template sync<DeviceType>();
  k_cost.template sync<DeviceType>();
  k_doExpansion.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCosineShiftExpCompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCosineShiftExpCompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCosineShiftExpCompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCosineShiftExpCompute<0,0> >(0,nanglelist),*this);
    }
  }

  if (eflag_global) energy += ev.evdwl;
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
void AngleCosineShiftExpKokkos<DeviceType>::operator()(TagAngleCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  const KK_FLOAT x20 = x(i2,0);
  const KK_FLOAT x21 = x(i2,1);
  const KK_FLOAT x22 = x(i2,2);

  // 1st bond

  const KK_FLOAT delx1 = x(i1,0) - x20;
  const KK_FLOAT dely1 = x(i1,1) - x21;
  const KK_FLOAT delz1 = x(i1,2) - x22;

  const KK_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const KK_FLOAT r1 = sqrt(rsq1);

  // 2nd bond

  const KK_FLOAT delx2 = x(i3,0) - x20;
  const KK_FLOAT dely2 = x(i3,1) - x21;
  const KK_FLOAT delz2 = x(i3,2) - x22;

  const KK_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const KK_FLOAT r2 = sqrt(rsq2);

  // c = cosine of angle, s = sine

  KK_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(1e-12)) s = static_cast<KK_FLOAT>(1e-12);

  const KK_FLOAT cccpsss = c*d_cost[type] + s*d_sint[type];
  const KK_FLOAT cssmscc = c*d_sint[type] - s*d_cost[type];

  const KK_FLOAT aa = d_a[type];
  const KK_FLOAT uumin = d_umin[type];

  KK_FLOAT eangle = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT ff;

  if (d_doExpansion[type]) {
    // Taylor expansion for small |a|
    if (eflag) eangle = -static_cast<KK_FLOAT>(0.125)*(static_cast<KK_FLOAT>(1.0)+cccpsss)
                        *(static_cast<KK_FLOAT>(4.0)+aa*(cccpsss-static_cast<KK_FLOAT>(1.0)))*uumin;
    ff = static_cast<KK_FLOAT>(0.25)*uumin*cssmscc*(static_cast<KK_FLOAT>(2.0)+aa*cccpsss)/s;
  } else {
    // full exponential formula
    const KK_FLOAT exp2 = exp(static_cast<KK_FLOAT>(0.5)*aa*(static_cast<KK_FLOAT>(1.0)+cccpsss));
    if (eflag) eangle = d_opt1[type]*(static_cast<KK_FLOAT>(1.0)-exp2);
    ff = static_cast<KK_FLOAT>(0.5)*aa*d_opt1[type]*exp2*cssmscc/s;
  }

  const KK_FLOAT a11 =  ff*c/rsq1;
  const KK_FLOAT a12 = -ff/(r1*r2);
  const KK_FLOAT a22 =  ff*c/rsq2;

  KK_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  // apply force to each of 3 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += f1[0];
    a_f(i1,1) += f1[1];
    a_f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= f1[0] + f3[0];
    a_f(i2,1) -= f1[1] + f3[1];
    a_f(i2,2) -= f1[2] + f3[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += f3[0];
    a_f(i3,1) += f3[1];
    a_f(i3,2) += f3[2];
  }

  if (EVFLAG) ev_tally(ev,i1,i2,i3,eangle,f1,f3,
                       delx1,dely1,delz1,delx2,dely2,delz2);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleCosineShiftExpKokkos<DeviceType>::operator()(TagAngleCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleCosineShiftExpCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCosineShiftExpKokkos<DeviceType>::allocate()
{
  AngleCosineShiftExp::allocate();

  int n = atom->nangletypes;
  k_umin        = DAT::tdual_kkfloat_1d("AngleCosineShiftExp::umin",n+1);
  k_a           = DAT::tdual_kkfloat_1d("AngleCosineShiftExp::a",n+1);
  k_opt1        = DAT::tdual_kkfloat_1d("AngleCosineShiftExp::opt1",n+1);
  k_sint        = DAT::tdual_kkfloat_1d("AngleCosineShiftExp::sint",n+1);
  k_cost        = DAT::tdual_kkfloat_1d("AngleCosineShiftExp::cost",n+1);
  k_doExpansion = DAT::tdual_int_1d("AngleCosineShiftExp::doExpansion",n+1);

  d_umin        = k_umin.template view<DeviceType>();
  d_a           = k_a.template view<DeviceType>();
  d_opt1        = k_opt1.template view<DeviceType>();
  d_sint        = k_sint.template view<DeviceType>();
  d_cost        = k_cost.template view<DeviceType>();
  d_doExpansion = k_doExpansion.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCosineShiftExpKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleCosineShiftExp::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_umin.view_host()[i]        = static_cast<KK_FLOAT>(umin[i]);
    k_a.view_host()[i]           = static_cast<KK_FLOAT>(a[i]);
    k_opt1.view_host()[i]        = static_cast<KK_FLOAT>(opt1[i]);
    k_sint.view_host()[i]        = static_cast<KK_FLOAT>(sint[i]);
    k_cost.view_host()[i]        = static_cast<KK_FLOAT>(cost[i]);
    k_doExpansion.view_host()[i] = doExpansion[i] ? 1 : 0;
  }

  k_umin.modify_host();
  k_a.modify_host();
  k_opt1.modify_host();
  k_sint.modify_host();
  k_cost.modify_host();
  k_doExpansion.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCosineShiftExpKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleCosineShiftExp::read_restart(fp);

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_umin.view_host()[i]        = static_cast<KK_FLOAT>(umin[i]);
    k_a.view_host()[i]           = static_cast<KK_FLOAT>(a[i]);
    k_opt1.view_host()[i]        = static_cast<KK_FLOAT>(opt1[i]);
    k_sint.view_host()[i]        = static_cast<KK_FLOAT>(sint[i]);
    k_cost.view_host()[i]        = static_cast<KK_FLOAT>(cost[i]);
    k_doExpansion.view_host()[i] = doExpansion[i] ? 1 : 0;
  }

  k_umin.modify_host();
  k_a.modify_host();
  k_opt1.modify_host();
  k_sint.modify_host();
  k_cost.modify_host();
  k_doExpansion.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleCosineShiftExpKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  KK_FLOAT eanglethird;
  KK_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += eangle;
      else {
        eanglethird = THIRD*eangle;

        if (i < nlocal) ev.evdwl += eanglethird;
        if (j < nlocal) ev.evdwl += eanglethird;
        if (k < nlocal) ev.evdwl += eanglethird;
      }
    }
    if (eflag_atom) {
      eanglethird = THIRD*eangle;

      if (newton_bond || i < nlocal) v_eatom[i] += eanglethird;
      if (newton_bond || j < nlocal) v_eatom[j] += eanglethird;
      if (newton_bond || k < nlocal) v_eatom[k] += eanglethird;
    }
  }

  if (vflag_either) {
    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (k < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += THIRD*v[0];
        v_vatom(i,1) += THIRD*v[1];
        v_vatom(i,2) += THIRD*v[2];
        v_vatom(i,3) += THIRD*v[3];
        v_vatom(i,4) += THIRD*v[4];
        v_vatom(i,5) += THIRD*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += THIRD*v[0];
        v_vatom(j,1) += THIRD*v[1];
        v_vatom(j,2) += THIRD*v[2];
        v_vatom(j,3) += THIRD*v[3];
        v_vatom(j,4) += THIRD*v[4];
        v_vatom(j,5) += THIRD*v[5];
      }
      if (newton_bond || k < nlocal) {
        v_vatom(k,0) += THIRD*v[0];
        v_vatom(k,1) += THIRD*v[1];
        v_vatom(k,2) += THIRD*v[2];
        v_vatom(k,3) += THIRD*v[3];
        v_vatom(k,4) += THIRD*v[4];
        v_vatom(k,5) += THIRD*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleCosineShiftExpKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleCosineShiftExpKokkos<LMPHostType>;
#endif
}
