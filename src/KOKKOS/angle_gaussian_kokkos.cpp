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

#include "angle_gaussian_kokkos.h"

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

static constexpr double SMALL  = 0.001;
static constexpr double SMALLG = 2.3e-308;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleGaussianKokkos<DeviceType>::AngleGaussianKokkos(LAMMPS *lmp) : AngleGaussian(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
  allocated_kokkos = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleGaussianKokkos<DeviceType>::~AngleGaussianKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleGaussianKokkos<DeviceType>::allocate_kokkos()
{
  int n = atom->nangletypes;

  if (!allocated_kokkos) {
    k_nterms = DAT::tdual_int_1d("AngleGaussian::nterms",n+1);
    k_alpha  = DAT::tdual_kkfloat_2d("AngleGaussian::alpha",n+1,nterms_max);
    k_width  = DAT::tdual_kkfloat_2d("AngleGaussian::width",n+1,nterms_max);
    k_theta0 = DAT::tdual_kkfloat_2d("AngleGaussian::theta0",n+1,nterms_max);
  } else {
    k_nterms.resize(n+1);
    k_alpha.resize(n+1,nterms_max);
    k_width.resize(n+1,nterms_max);
    k_theta0.resize(n+1,nterms_max);
  }

  d_nterms = k_nterms.template view<DeviceType>();
  d_alpha  = k_alpha.template view<DeviceType>();
  d_width  = k_width.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();

  allocated_kokkos = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleGaussianKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

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

  k_nterms.template sync<DeviceType>();
  k_angle_temperature.template sync<DeviceType>();
  k_alpha.template sync<DeviceType>();
  k_width.template sync<DeviceType>();
  k_theta0.template sync<DeviceType>();

  boltz = static_cast<KK_FLOAT>(force->boltz);

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleGaussianCompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleGaussianCompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleGaussianCompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleGaussianCompute<0,0> >(0,nanglelist),*this);
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
void AngleGaussianKokkos<DeviceType>::operator()(TagAngleGaussianCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  const KK_FLOAT delx1 = x(i1,0) - x(i2,0);
  const KK_FLOAT dely1 = x(i1,1) - x(i2,1);
  const KK_FLOAT delz1 = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const KK_FLOAT r1 = sqrt(rsq1);

  const KK_FLOAT delx2 = x(i3,0) - x(i2,0);
  const KK_FLOAT dely2 = x(i3,1) - x(i2,1);
  const KK_FLOAT delz2 = x(i3,2) - x(i2,2);

  const KK_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const KK_FLOAT r2 = sqrt(rsq2);

  KK_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  s = static_cast<KK_FLOAT>(1.0)/s;

  const KK_FLOAT theta = acos(c);

  KK_FLOAT sum_g_i = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT sum_numerator = static_cast<KK_FLOAT>(0.0);
  const int nt = d_nterms[type];
  for (int i = 0; i < nt; i++) {
    const KK_FLOAT dtheta    = theta - d_theta0(type,i);
    const KK_FLOAT w         = d_width(type,i);
    const KK_FLOAT prefactor = d_alpha(type,i) / (w * sqrt(static_cast<KK_FLOAT>(MY_PI2)));
    const KK_FLOAT exponent  = static_cast<KK_FLOAT>(-2.0) * dtheta * dtheta / (w * w);
    const KK_FLOAT g_i       = prefactor * exp(exponent);
    sum_g_i       += g_i;
    sum_numerator += g_i * dtheta / (w * w);
  }

  // avoid overflow
  if (sum_g_i < sum_numerator * static_cast<KK_FLOAT>(SMALLG))
    sum_g_i = sum_numerator * static_cast<KK_FLOAT>(SMALLG);

  const KK_FLOAT kbt = boltz * d_angle_temperature[type];

  KK_FLOAT eangle = static_cast<KK_FLOAT>(0.0);
  if (eflag) eangle = -kbt * log(sum_g_i);

  const KK_FLOAT a   = static_cast<KK_FLOAT>(-4.0) * kbt * (sum_numerator / sum_g_i) * s;
  const KK_FLOAT a11 = a*c / rsq1;
  const KK_FLOAT a12 = -a / (r1*r2);
  const KK_FLOAT a22 = a*c / rsq2;

  KK_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += static_cast<KK_ACC_FLOAT>(f1[0]);
    a_f(i1,1) += static_cast<KK_ACC_FLOAT>(f1[1]);
    a_f(i1,2) += static_cast<KK_ACC_FLOAT>(f1[2]);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= static_cast<KK_ACC_FLOAT>(f1[0] + f3[0]);
    a_f(i2,1) -= static_cast<KK_ACC_FLOAT>(f1[1] + f3[1]);
    a_f(i2,2) -= static_cast<KK_ACC_FLOAT>(f1[2] + f3[2]);
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += static_cast<KK_ACC_FLOAT>(f3[0]);
    a_f(i3,1) += static_cast<KK_ACC_FLOAT>(f3[1]);
    a_f(i3,2) += static_cast<KK_ACC_FLOAT>(f3[2]);
  }

  if (EVFLAG) ev_tally(ev,i1,i2,i3,eangle,f1,f3,
                       delx1,dely1,delz1,delx2,dely2,delz2);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleGaussianKokkos<DeviceType>::operator()(TagAngleGaussianCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleGaussianCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleGaussianKokkos<DeviceType>::allocate()
{
  AngleGaussian::allocate();

  int n = atom->nangletypes;
  k_nterms            = DAT::tdual_int_1d("AngleGaussian::nterms",n+1);
  k_angle_temperature = DAT::tdual_kkfloat_1d("AngleGaussian::angle_temperature",n+1);

  d_nterms            = k_nterms.template view<DeviceType>();
  d_angle_temperature = k_angle_temperature.template view<DeviceType>();
  // 2D views are created lazily in allocate_kokkos() once max_nterms is known
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleGaussianKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleGaussian::coeff(narg, arg);
  allocate_kokkos();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->ndihedraltypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_nterms.view_host()[i]             = nterms[i];
    k_angle_temperature.view_host()[i]  = angle_temperature[i];
    for (int j = 0; j < nterms[i]; j++) {
      k_alpha.view_host()(i,j)  = alpha[i][j];
      k_width.view_host()(i,j)  = width[i][j];
      k_theta0.view_host()(i,j) = theta0[i][j];
    }
  }

  k_nterms.modify_host();
  k_angle_temperature.modify_host();
  k_alpha.modify_host();
  k_width.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleGaussianKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleGaussian::read_restart(fp);
  allocate_kokkos();

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_nterms.view_host()[i]             = nterms[i];
    k_angle_temperature.view_host()[i]  = angle_temperature[i];
    for (int j = 0; j < nterms[i]; j++) {
      k_alpha.view_host()(i,j)  = alpha[i][j];
      k_width.view_host()(i,j)  = width[i][j];
      k_theta0.view_host()(i,j) = theta0[i][j];
    }
  }

  k_nterms.modify_host();
  k_angle_temperature.modify_host();
  k_alpha.modify_host();
  k_width.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleGaussianKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(eangle);
      else {
        KK_ACC_FLOAT eanglethird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*eangle);
        if (i < nlocal) ev.evdwl += eanglethird;
        if (j < nlocal) ev.evdwl += eanglethird;
        if (k < nlocal) ev.evdwl += eanglethird;
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT eanglethird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*eangle);
      if (newton_bond || i < nlocal) v_eatom[i] += static_cast<KK_ACC_FLOAT>(eanglethird);
      if (newton_bond || j < nlocal) v_eatom[j] += static_cast<KK_ACC_FLOAT>(eanglethird);
      if (newton_bond || k < nlocal) v_eatom[k] += static_cast<KK_ACC_FLOAT>(eanglethird);
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v_third_acc[6];
    v_third_acc[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[0] + delx2*f3[0]));
    v_third_acc[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(dely1*f1[1] + dely2*f3[1]));
    v_third_acc[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delz1*f1[2] + delz2*f3[2]));
    v_third_acc[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[1] + delx2*f3[1]));
    v_third_acc[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[2] + delx2*f3[2]));
    v_third_acc[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(dely1*f1[2] + dely2*f3[2]));

    if (vflag_global) {
      if (newton_bond) {
        for (int n = 0; n < 6; n++)
          ev.v[n] += static_cast<KK_ACC_FLOAT>(static_cast<KK_ACC_FLOAT>(3.0)*v_third_acc[n]);
      } else {
        if (i < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
        if (j < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
        if (k < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(i,n) += v_third_acc[n];
      }
      if (newton_bond || j < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(j,n) += v_third_acc[n];
      }
      if (newton_bond || k < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(k,n) += v_third_acc[n];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleGaussianKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleGaussianKokkos<LMPHostType>;
#endif
}
