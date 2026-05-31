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
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

// pre-kernel: compute the M (virtual charge) site for every O atom
struct TagPairLJCutTIP4PCutNewsite{};

// main kernel: templated on energy/virial flag
template<int EVFLAG>
struct TagPairLJCutTIP4PCutCompute{};

template<class DeviceType>
class PairLJCutTIP4PCutKokkos : public PairLJCutTIP4PCut {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairLJCutTIP4PCutKokkos(class LAMMPS *);
  ~PairLJCutTIP4PCutKokkos() override;

  void compute(int, int) override;
  void init_style() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PCutNewsite, const int&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PCutCompute<EVFLAG>, const int&, EV_FLOAT&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PCutCompute<EVFLAG>, const int&) const;

 protected:
  // find the periodic image of j closest to i (walks the sametag chain)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int closest_image(const int, int) const;

  // standard pairwise (LJ) energy/virial tally for a half neighbor list
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &evdwl,
                const KK_FLOAT &fpair, const KK_FLOAT &delx, const KK_FLOAT &dely,
                const KK_FLOAT &delz) const;

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

  // per-type-pair LJ coefficients (device)
  typename AT::t_kkfloat_2d d_lj1, d_lj2, d_lj3, d_lj4, d_offset, d_cut_ljsq;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  int newton_pair, neighflag;
  int nlocal, nall, eflag, vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT m_alpha, m_qdist;
  KK_FLOAT m_cut_coulsq, m_cut_coulsqplus;
  int m_typeO, m_typeH;
};

}    // namespace LAMMPS_NS

#endif
#endif
