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
PairStyle(gran/hooke/history/kk,PairGranHookeHistoryKokkos<LMPDeviceType>);
PairStyle(gran/hooke/history/kk/device,PairGranHookeHistoryKokkos<LMPDeviceType>);
PairStyle(gran/hooke/history/kk/host,PairGranHookeHistoryKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_KOKKOS_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_KOKKOS_H

#include "pair_gran_hooke_history.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixNeighHistoryKokkos;

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int SHEARUPDATE>
struct TagPairGranHookeHistoryCompute {};

template <class DeviceType>
class PairGranHookeHistoryKokkos : public PairGranHookeHistory {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairGranHookeHistoryKokkos(class LAMMPS *);
  ~PairGranHookeHistoryKokkos() override;
  void compute(int, int) override;
  void init_style() override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int SHEARUPDATE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int SHEARUPDATE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,SHEARUPDATE>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
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
  typename AT::t_kkfloat_1d_3 torque;
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

  typename AT::t_int_2d d_firsttouch;
  typename AT::t_kkfloat_2d d_firstshear;

  typename AT::t_neighbors_2d d_neighbors_touch;
  typename AT::t_int_1d d_numneigh_touch;

  int newton_pair;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const {return j >> SBBITS & 3;}

  friend void pair_virial_fdotr_compute<PairGranHookeHistoryKokkos>(PairGranHookeHistoryKokkos*);
};

}

#endif
#endif

