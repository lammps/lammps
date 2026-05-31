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
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

// pre-kernel: compute the M (virtual charge) site for every O atom
struct TagPairTIP4PLongNewsite{};

// main kernel: templated on energy/virial flag
template<int EVFLAG>
struct TagPairTIP4PLongCompute{};

template<class DeviceType>
class PairTIP4PLongKokkos : public PairTIP4PLong {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairTIP4PLongKokkos(class LAMMPS *);
  ~PairTIP4PLongKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  void init_tables(double cut_coul, double *cut_respa) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PLongNewsite, const int&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PLongCompute<EVFLAG>, const int&, EV_FLOAT&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTIP4PLongCompute<EVFLAG>, const int&) const;

 protected:
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int closest_image(const int, int) const;

  // M-site (TIP4P) energy/virial tally; key encodes whether i and/or j are O
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_tip4p(EV_FLOAT &ev, const int &key, const int (&vlist)[6],
                      const KK_FLOAT (&v)[6], const KK_FLOAT &ecoul) const;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d_randomread tag;
  typename AT::t_int_1d d_sametag;

  typename AT::t_kkfloat_1d_3 d_newsite;
  typename AT::t_int_1d_3 d_hneigh;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  // Ewald real-space + optional coulomb interpolation tables
  typename AT::t_kkfloat_1d d_rtable, d_drtable, d_ftable, d_dftable,
                            d_ctable, d_dctable, d_etable, d_detable;
  KK_FLOAT g_ewald_kk, tabinnersq_kk;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  int newton_pair, neighflag;
  int nlocal, nall, eflag, vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT m_alpha, m_qdist;
  KK_FLOAT m_cut_coulsq, m_cut_coulsqplus;
  int m_typeO, m_typeH;
  int m_ncoultablebits, m_ncoulmask, m_ncoulshiftbits;
};

}    // namespace LAMMPS_NS

#endif
#endif
