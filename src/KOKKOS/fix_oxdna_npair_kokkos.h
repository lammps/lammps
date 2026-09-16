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

#ifdef FIX_CLASS
// clang-format off
FixStyle(OXDNA/NPAIR/kk,FixOxdnaNpairKokkos<LMPDeviceType>);
FixStyle(OXDNA/NPAIR/kk/device,FixOxdnaNpairKokkos<LMPDeviceType>);
FixStyle(OXDNA/NPAIR/kk/host,FixOxdnaNpairKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_OXDNA_NPAIR_KOKKOS_H
#define LMP_FIX_OXDNA_NPAIR_KOKKOS_H

#include "fix.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixOxdnaNpairNeighScreen{};
struct TagFixOxdnaNpairFill{};

template<class DeviceType>
class FixOxdnaNpairKokkos : public Fix {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixOxdnaNpairKokkos(class LAMMPS *, int, char **);
  ~FixOxdnaNpairKokkos() override;

  void init() override;
  int setmask() override;
  void init_list(int, class NeighList *) override;
  void min_setup_pre_force(int);
  void min_pre_force(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  void compute_neigh_screen_to_npair();

  // Derived COM screen cutoff. Each consuming pair style (hbond / xstk /
  // coaxstk) registers its max interaction-site cutoff plus the site-offset
  // margin (so a center-of-mass test never drops an interacting pair) in its
  // init_one; the screen uses the largest request. Falls back to the historical
  // r < 2.0 if nothing registers.
  void request_screen_cutoff(double cut_com) {
    if (cut_com > screen_cut_max) screen_cut_max = cut_com;
  }

  // By default, screened-pair rebuild is device-only. This override lets
  // selected styles force npair rebuilds on host backends too.
  // Right now, only the oxdna3/xstk/kk pair style makes use of this - this was a design choice
  // due to the need for pre-atom_mapping in KOKKOS which relies on a given neighbor list.
  // Seems a bit daft to calculate this for both a 'normal' neighbor list
  // and our custom oxdna npair list which would carry unwanted overhead.
  // Hence, the option to force the rebuild on all backends.
  void set_force_screening_all_backends(bool value) { force_screening_all_backends = value; }

  // Direct packed (a, b) pair lookup for coalesced access on GPUs.
  DAT::tdual_uint64_1d k_pairs_screened;
  typename AT::t_uint64_1d d_pairs_screened;
  int screened_pair_count; // ComputeGPUPair functors use this in place of the usual anum.

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaNpairNeighScreen, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaNpairFill, const int &) const;

 private:

  class NeighList *list;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;

  int anum;
  int neighflag, last_allocate;
  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;
  // Screening takes place on GPUs only
  DAT::tdual_int_1d k_numneigh_screened;
  typename AT::t_int_1d d_numneigh_screened;
  DAT::tdual_int_1d k_screened_offsets;
  typename AT::t_int_1d d_screened_offsets;
  DAT::tdual_int_scalar k_screened_pair_count;
  typename AT::t_int_scalar d_screened_pair_count;
  int screened_max_atoms;
  int screened_max_neigh;
  double screen_cut_max;   // max COM screen cutoff requested by consuming styles (host)
  KK_FLOAT screen_cutsq;   // screen_cut_max^2, read on device by screen_pair_fast
  bool force_screening_all_backends;

  void update_screen_cutsq();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool screen_pair_fast(const int &braw,
    const KK_FLOAT &a_com0, const KK_FLOAT &a_com1, const KK_FLOAT &a_com2) const;
};

}    // namespace LAMMPS_NS
#endif
#endif
