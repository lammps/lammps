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
PairStyle(lj/cut/tip4p/cut/kk,PairLJCutTIP4PCutKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/cut/kk/device,PairLJCutTIP4PCutKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/cut/kk/host,PairLJCutTIP4PCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_TIP4P_CUT_KOKKOS_H
#define LMP_PAIR_LJ_CUT_TIP4P_CUT_KOKKOS_H

#include "pair_lj_cut_tip4p_cut.h"
#include "pair_tip4p_kokkos.h"

namespace LAMMPS_NS {

template<int EVFLAG>
struct TagPairLJCutTIP4PCutCompute{};

template<class DeviceType>
class PairLJCutTIP4PCutKokkos : public PairTIP4PKokkos<DeviceType,PairLJCutTIP4PCut> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typedef PairTIP4PKokkos<DeviceType,PairLJCutTIP4PCut> Base;

  PairLJCutTIP4PCutKokkos(class LAMMPS *lmp) : Base(lmp) {}

  void compute(int, int) override;
  double init_one(int, int) override;

  using Base::operator();

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PCutCompute<EVFLAG>, const int&, EV_FLOAT&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PCutCompute<EVFLAG>, const int&) const;

 protected:
  // per-type-pair LJ coefficients (device); re-synced when init_one() ran
  typename AT::t_kkfloat_2d d_lj1, d_lj2, d_lj3, d_lj4, d_offset, d_cut_ljsq;
  int lj_coeffs_synced = 0;

  using Base::x; using Base::f; using Base::q; using Base::type;
  using Base::d_newsite; using Base::d_hneigh; using Base::d_h_missing;
  using Base::d_neighbors; using Base::d_numneigh; using Base::d_ilist;
  using Base::m_typeO;
  using Base::m_cut_coulsq; using Base::m_cut_coulsqplus;
  using Base::qqrd2e; using Base::special_coul; using Base::special_lj;
  using Base::sbmask; using Base::ev_tally_tip4p; using Base::apply_site_force;
  using Base::ev_tally;
};

}    // namespace LAMMPS_NS

#endif
#endif
