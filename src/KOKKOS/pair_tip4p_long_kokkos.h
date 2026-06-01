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
PairStyle(tip4p/long/kk,PairTIP4PLongKokkos<LMPDeviceType>);
PairStyle(tip4p/long/kk/device,PairTIP4PLongKokkos<LMPDeviceType>);
PairStyle(tip4p/long/kk/host,PairTIP4PLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TIP4P_LONG_KOKKOS_H
#define LMP_PAIR_TIP4P_LONG_KOKKOS_H

#include "pair_tip4p_long.h"
#include "pair_tip4p_kokkos.h"

namespace LAMMPS_NS {

template<int EVFLAG>
struct TagPairTIP4PLongCompute{};

template<class DeviceType>
class PairTIP4PLongKokkos : public PairTIP4PKokkos<DeviceType,PairTIP4PLong> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typedef PairTIP4PKokkos<DeviceType,PairTIP4PLong> Base;

  PairTIP4PLongKokkos(class LAMMPS *lmp) : Base(lmp) {}

  void compute(int, int) override;
  void init_tables(double cut_coul, double *cut_respa) override;

  using Base::operator();

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PLongCompute<EVFLAG>, const int&, EV_FLOAT&) const;

 protected:
  // Ewald real-space + optional coulomb interpolation tables (device)
  typename AT::t_kkfloat_1d d_rtable, d_drtable, d_ftable, d_dftable,
                            d_ctable, d_dctable, d_etable, d_detable;
  KK_FLOAT g_ewald_kk, tabinnersq_kk;
  int m_ncoultablebits, m_ncoulmask, m_ncoulshiftbits;

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
