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
AngleStyle(class2/p6/kk,AngleClass2P6Kokkos<LMPDeviceType>);
AngleStyle(class2/p6/kk/device,AngleClass2P6Kokkos<LMPDeviceType>);
AngleStyle(class2/p6/kk/host,AngleClass2P6Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_CLASS2_P6_KOKKOS_H
#define LMP_ANGLE_CLASS2_P6_KOKKOS_H

#include "angle_class2_p6.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleClass2P6Compute{};

template<class DeviceType>
class AngleClass2P6Kokkos : public AngleClass2P6 {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  AngleClass2P6Kokkos(class LAMMPS *);
  ~AngleClass2P6Kokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2P6Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2P6Compute<NEWTON_BOND,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_2d_lr anglelist;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_kkfloat_1d k_theta0;
  DAT::tdual_kkfloat_1d k_k2, k_k3, k_k4, k_k5, k_k6;
  DAT::tdual_kkfloat_1d k_bb_k, k_bb_r1, k_bb_r2;
  DAT::tdual_kkfloat_1d k_ba_k1, k_ba_k2, k_ba_r1, k_ba_r2;

  typename AT::t_kkfloat_1d d_theta0;
  typename AT::t_kkfloat_1d d_k2, d_k3, d_k4, d_k5, d_k6;
  typename AT::t_kkfloat_1d d_bb_k, d_bb_r1, d_bb_r2;
  typename AT::t_kkfloat_1d d_ba_k1, d_ba_k2, d_ba_r1, d_ba_r2;

  void allocate() override;
};

}

#endif
#endif
