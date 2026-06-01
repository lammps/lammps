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
PairStyle(tip4p/cut/kk,PairTIP4PCutKokkos<LMPDeviceType>);
PairStyle(tip4p/cut/kk/device,PairTIP4PCutKokkos<LMPDeviceType>);
PairStyle(tip4p/cut/kk/host,PairTIP4PCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TIP4P_CUT_KOKKOS_H
#define LMP_PAIR_TIP4P_CUT_KOKKOS_H

#include "pair_tip4p_cut.h"
#include "pair_tip4p_kokkos.h"

namespace LAMMPS_NS {

template<int EVFLAG>
struct TagPairTIP4PCutCompute{};

template<class DeviceType>
class PairTIP4PCutKokkos : public PairTIP4PKokkos<DeviceType,PairTIP4PCut> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typedef PairTIP4PKokkos<DeviceType,PairTIP4PCut> Base;

  PairTIP4PCutKokkos(class LAMMPS *lmp) : Base(lmp) {}

  void compute(int, int) override;

  // un-hide the base-class M-site pre-kernel operator()
  using Base::operator();

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PCutCompute<EVFLAG>, const int&, EV_FLOAT&) const;

 protected:
  using Base::x; using Base::f; using Base::q; using Base::type;
  using Base::d_newsite; using Base::d_hneigh;
  using Base::d_neighbors; using Base::d_numneigh; using Base::d_ilist;
  using Base::m_typeO; using Base::m_alpha;
  using Base::m_cut_coulsq; using Base::m_cut_coulsqplus;
  using Base::qqrd2e; using Base::special_coul;
  using Base::sbmask; using Base::ev_tally_tip4p;
};

}    // namespace LAMMPS_NS

#endif
#endif
