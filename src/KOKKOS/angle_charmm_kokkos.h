/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(charmm/kk,AngleCharmmKokkos<LMPDeviceType>);
AngleStyle(charmm/kk/device,AngleCharmmKokkos<LMPDeviceType>);
AngleStyle(charmm/kk/host,AngleCharmmKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_CHARMM_KOKKOS_H
#define LMP_ANGLE_CHARMM_KOKKOS_H

#include "angle_charmm.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleCharmmCompute{};

template<class DeviceType>
class AngleCharmmKokkos : public AngleCharmm {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  AngleCharmmKokkos(class LAMMPS *);
  ~AngleCharmmKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const;

  using KKDeviceType = typename KKDevice<DeviceType>::value;
  TransformView<KK_ACC_FLOAT*,double*,Kokkos::LayoutRight,KKDeviceType> k_eatom;
  TransformView<KK_ACC_FLOAT*[6],double*[6],LMPDeviceLayout,KKDeviceType> k_vatom;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename Kokkos::View<KK_ACC_FLOAT*[3],DAT::t_kkacc_1d_3::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > f;
  typename AT::t_int_2d_lr anglelist;
  Kokkos::View<KK_ACC_FLOAT*,Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic>> d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6],LMPDeviceLayout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic>> d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  typename AT::t_kkfloat_1d d_k;
  typename AT::t_kkfloat_1d d_theta0;
  typename AT::t_kkfloat_1d d_k_ub;
  typename AT::t_kkfloat_1d d_r_ub;

  void allocate() override;
};

}

#endif
#endif

