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
AngleStyle(gaussian/kk,AngleGaussianKokkos<LMPDeviceType>);
AngleStyle(gaussian/kk/device,AngleGaussianKokkos<LMPDeviceType>);
AngleStyle(gaussian/kk/host,AngleGaussianKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_GAUSSIAN_KOKKOS_H
#define LMP_ANGLE_GAUSSIAN_KOKKOS_H

#include "angle_gaussian.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleGaussianCompute{};

template<class DeviceType>
class AngleGaussianKokkos : public AngleGaussian {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  AngleGaussianKokkos(class LAMMPS *);
  ~AngleGaussianKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleGaussianCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleGaussianCompute<NEWTON_BOND,EVFLAG>, const int&) const;

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
  KK_FLOAT boltz;

  DAT::tdual_int_1d k_nterms;
  DAT::tdual_kkfloat_1d k_angle_temperature;
  DAT::tdual_kkfloat_2d k_alpha;
  DAT::tdual_kkfloat_2d k_width;
  DAT::tdual_kkfloat_2d k_theta0;

  typename AT::t_int_1d d_nterms;
  typename AT::t_kkfloat_1d d_angle_temperature;
  typename AT::t_kkfloat_2d d_alpha;
  typename AT::t_kkfloat_2d d_width;
  typename AT::t_kkfloat_2d d_theta0;

  void allocate() override;
  void allocate_kokkos();
  int allocated_kokkos;
};
}
#endif
#endif
