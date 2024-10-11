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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(brownian/kk,PairBrownianKokkos<LMPDeviceType>);
PairStyle(brownian/kk/device,PairBrownianKokkos<LMPDeviceType>);
PairStyle(brownian/kk/host,PairBrownianKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BROWNIAN_KOKKOS_H
#define LMP_PAIR_BROWNIAN_KOKKOS_H

#include "pair_brownian.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "Kokkos_Random.hpp"
#include "comm_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
struct TagPairBrownianCompute {};

template<class DeviceType>
class PairBrownianKokkos : public PairBrownian, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairBrownianKokkos(class LAMMPS *);
  ~PairBrownianKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBrownianCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBrownianCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                    F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                    X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;

 protected:
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread radius;

  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int newton_pair;
  double special_lj[4];

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;
  typename AT::tdual_ffloat_2d k_cut_inner;
  typename AT::t_ffloat_2d d_cut_inner;

  int neighflag;
  int nlocal,nall,eflag,vflag;
  LMP_FLOAT vxmu2f;

  LMP_FLOAT prethermostat;

  void allocate() override;

  KOKKOS_INLINE_FUNCTION
  void set_3_orthogonal_vectors(const double p1[3], double * const p2, double * const p3) const {
    double norm;
    int ix, iy, iz;

    // find the index of maximum magnitude and store it in iz

    if (fabs(p1[0]) > fabs(p1[1])) {
      iz = 0;
      ix = 1;
      iy = 2;
    } else {
      iz = 1;
      ix = 2;
      iy = 0;
    }

    if (iz == 0) {
      if (fabs(p1[0]) < fabs(p1[2])) {
        iz = 2;
        ix = 0;
        iy = 1;
      }
    } else {
      if (fabs(p1[1]) < fabs(p1[2])) {
        iz = 2;
        ix = 0;
        iy = 1;
      }
    }

    // set p2 arbitrarily such that it's orthogonal to p1

    p2[ix] = 1.0;
    p2[iy] = 1.0;
    p2[iz] = -(p1[ix] * p2[ix] + p1[iy] * p2[iy]) / p1[iz];

    // normalize p2

    norm = sqrt(p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);

    p2[0] = p2[0] / norm;
    p2[1] = p2[1] / norm;
    p2[2] = p2[2] / norm;

    // Set p3 by taking the cross product p3=p2xp1

    p3[0] = p1[1] * p2[2] - p1[2] * p2[1];
    p3[1] = p1[2] * p2[0] - p1[0] * p2[2];
    p3[2] = p1[0] * p2[1] - p1[1] * p2[0];
  };

  friend void pair_virial_fdotr_compute<PairBrownianKokkos>(PairBrownianKokkos*);

  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
};

}

#endif
#endif

