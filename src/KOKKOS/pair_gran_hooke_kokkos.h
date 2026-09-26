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
PairStyle(gran/hooke/kk,PairGranHookeKokkos<LMPDeviceType>);
PairStyle(gran/hooke/kk/device,PairGranHookeKokkos<LMPDeviceType>);
PairStyle(gran/hooke/kk/host,PairGranHookeKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRAN_HOOKE_KOKKOS_H
#define LMP_PAIR_GRAN_HOOKE_KOKKOS_H

#include "pair_gran_hooke.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG>
struct TagPairGranHookeCompute {};

template <class DeviceType>
class PairGranHookeKokkos : public PairGranHooke {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairGranHookeKokkos(class LAMMPS *);
  ~PairGranHookeKokkos() override;
  void compute(int, int) override;
  void init_style() override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkfloat_1d_3_randomread omega;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread radius;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;



  int newton_pair;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT kt_kk, kn_kk, xmu_kk, gammat_kk, gamman_kk, dt_kk;


// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const {return j >> SBBITS & 3;}

  friend void pair_virial_fdotr_compute<PairGranHookeKokkos>(PairGranHookeKokkos*);
};

}

#endif
#endif

