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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "angle_table_kokkos.h"

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

static constexpr double SMALL = 0.001;
static constexpr int LINEAR_STYLE = 0;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleTableKokkos<DeviceType>::AngleTableKokkos(LAMMPS *lmp) : AngleTable(lmp)
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
AngleTableKokkos<DeviceType>::~AngleTableKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  //atomKK->sync(execution_space,datamask_read);
  //  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  //  else atomKK->modified(execution_space,F_MASK);

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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleTableCompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleTableCompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleTableCompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleTableCompute<0,0> >(0,nanglelist),*this);
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
void AngleTableKokkos<DeviceType>::operator()(TagAngleTableCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  // 1st bond

  const KK_FLOAT delx1 = x(i1,0) - x(i2,0);
  const KK_FLOAT dely1 = x(i1,1) - x(i2,1);
  const KK_FLOAT delz1 = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const KK_FLOAT r1 = Kokkos::sqrt(rsq1);

  // 2nd bond

  const KK_FLOAT delx2 = x(i3,0) - x(i2,0);
  const KK_FLOAT dely2 = x(i3,1) - x(i2,1);
  const KK_FLOAT delz2 = x(i3,2) - x(i2,2);

  const KK_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const KK_FLOAT r2 = Kokkos::sqrt(rsq2);

  // angle (cos and sin)

  KK_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = Kokkos::sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  s = static_cast<KK_FLOAT>(1.0)/s;

  // tabulated force & energy

  KK_FLOAT u,mdu;
  uf_lookup_kk(type,Kokkos::acos(c),u,mdu);

  KK_FLOAT eangle = 0;
  if (eflag) eangle = u;

  const KK_FLOAT a = mdu * s;
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

  // apply force to each of 3 atoms

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
void AngleTableKokkos<DeviceType>::operator()(TagAngleTableCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleTableCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ----------------------------------------------------------------------
   copy the tabulated data built on the host by compute_table() into
   device views.  this runs once per setup, after every coeff() call, so
   restarts and repeated angle_coeff commands are covered
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleTableKokkos<DeviceType>::setup_tables()
{
  const int n = atom->nangletypes;

  k_tabindex = DAT::tdual_int_1d("AngleTable::tabindex",n+1);
  k_invdelta = DAT::tdual_kkfloat_1d("AngleTable::invdelta",ntables);
  k_deltasq6 = DAT::tdual_kkfloat_1d("AngleTable::deltasq6",ntables);

  k_ang = DAT::tdual_kkfloat_2d("AngleTable::ang",ntables,tablength+1);
  k_e = DAT::tdual_kkfloat_2d("AngleTable::e",ntables,tablength+1);
  k_de = DAT::tdual_kkfloat_2d("AngleTable::de",ntables,tablength+1);
  k_f = DAT::tdual_kkfloat_2d("AngleTable::f",ntables,tablength+1);
  k_df = DAT::tdual_kkfloat_2d("AngleTable::df",ntables,tablength+1);
  k_e2 = DAT::tdual_kkfloat_2d("AngleTable::e2",ntables,tablength+1);
  k_f2 = DAT::tdual_kkfloat_2d("AngleTable::f2",ntables,tablength+1);

  for (int i = 1; i <= n; i++) k_tabindex.view_host()(i) = tabindex[i];

  for (int m = 0; m < ntables; m++) {
    const Table *tb = &tables[m];
    k_invdelta.view_host()(m) = static_cast<KK_FLOAT>(tb->invdelta);
    k_deltasq6.view_host()(m) = static_cast<KK_FLOAT>(tb->deltasq6);
    for (int i = 0; i < tablength; i++) {
      k_ang.view_host()(m,i) = static_cast<KK_FLOAT>(tb->ang[i]);
      k_e.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e[i]);
      k_de.view_host()(m,i) = static_cast<KK_FLOAT>(tb->de[i]);
      k_f.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f[i]);
      k_df.view_host()(m,i) = static_cast<KK_FLOAT>(tb->df[i]);
      k_e2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->e2[i]);
      k_f2.view_host()(m,i) = static_cast<KK_FLOAT>(tb->f2[i]);
    }
  }

  k_tabindex.modify_host(); k_tabindex.template sync<DeviceType>();
  k_invdelta.modify_host(); k_invdelta.template sync<DeviceType>();
  k_deltasq6.modify_host(); k_deltasq6.template sync<DeviceType>();
  k_ang.modify_host(); k_ang.template sync<DeviceType>();
  k_e.modify_host(); k_e.template sync<DeviceType>();
  k_de.modify_host(); k_de.template sync<DeviceType>();
  k_f.modify_host(); k_f.template sync<DeviceType>();
  k_df.modify_host(); k_df.template sync<DeviceType>();
  k_e2.modify_host(); k_e2.template sync<DeviceType>();
  k_f2.modify_host(); k_f2.template sync<DeviceType>();

  d_tabindex = k_tabindex.template view<DeviceType>();
  d_invdelta = k_invdelta.template view<DeviceType>();
  d_deltasq6 = k_deltasq6.template view<DeviceType>();
  d_ang = k_ang.template view<DeviceType>();
  d_e = k_e.template view<DeviceType>();
  d_de = k_de.template view<DeviceType>();
  d_f = k_f.template view<DeviceType>();
  d_df = k_df.template view<DeviceType>();
  d_e2 = k_e2.template view<DeviceType>();
  d_f2 = k_f2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleTableKokkos<DeviceType>::init_style()
{
  AngleTable::init_style();

  setup_tables();
}

/* ----------------------------------------------------------------------
   device version of AngleTable::uf_lookup()
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleTableKokkos<DeviceType>::uf_lookup_kk(const int &type, const KK_FLOAT &x_in,
                                                KK_FLOAT &u, KK_FLOAT &mdu) const
{
  const int tb = d_tabindex[type];
  const KK_FLOAT invdelta = d_invdelta[tb];

  // invdelta is based on tablength-1

  int itable = static_cast<int>(x_in * invdelta);
  if (itable < 0) itable = 0;
  if (itable >= tablength) itable = tablength - 1;

  const KK_FLOAT b = (x_in - d_ang(tb,itable)) * invdelta;

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

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleTableKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  // The eatom and vatom arrays are atomic
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
template class AngleTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleTableKokkos<LMPHostType>;
#endif
}

