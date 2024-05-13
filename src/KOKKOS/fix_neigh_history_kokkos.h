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
FixStyle(NEIGH_HISTORY/KK,FixNeighHistoryKokkos<LMPDeviceType>);
FixStyle(NEIGH_HISTORY/KK/DEVICE,FixNeighHistoryKokkos<LMPDeviceType>);
FixStyle(NEIGH_HISTORY/KK/HOST,FixNeighHistoryKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NEIGH_HISTORY_KOKKOS_H
#define LMP_FIX_NEIGH_HISTORY_KOKKOS_H

#include "fix_neigh_history.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

struct TagFixNeighHistoryPreExchange{};
struct TagFixNeighHistoryPostNeighbor{};
struct TagFixNeighHistoryPackExchange{};
struct TagFixNeighHistoryUnpackExchange{};

template <class DeviceType>
class FixNeighHistoryKokkos : public FixNeighHistory, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef int value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNeighHistoryKokkos(class LAMMPS *, int, char **);
  ~FixNeighHistoryKokkos() override;

  void pre_exchange() override;
  void post_neighbor() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  double memory_usage() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNeighHistoryPreExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNeighHistoryPostNeighbor, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNeighHistoryPackExchange, const int&, int &, const bool &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNeighHistoryUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
			   DAT::tdual_int_1d k_sendlist,
			   DAT::tdual_int_1d k_copylist,
			   ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nrecv1,int nrecv1extra,
                              ExecutionSpace space) override;

  typename DAT::tdual_int_2d k_firstflag;
  typename DAT::tdual_float_2d k_firstvalue;

 private:
  int nrecv1,nextrarecv1;
  int nlocal,nsend,beyond_contact;

  typename AT::t_tagint_1d tag;

  typename AT::t_int_2d d_firstflag;
  typename AT::t_float_2d d_firstvalue;

  DAT::tdual_int_1d k_npartner;
  DAT::tdual_tagint_2d k_partner;
  DAT::tdual_float_2d k_valuepartner;

  typename AT::t_int_1d d_npartner;
  typename AT::t_tagint_2d d_partner;
  typename AT::t_float_2d d_valuepartner;

  typename AT::t_int_1d d_sendlist;
  typename AT::t_xfloat_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_scalar d_resize,d_count;
  HAT::t_int_scalar h_resize,h_count;

  void pre_exchange_no_newton() override;

  // Shift by HISTBITS and check the first bit
  KOKKOS_INLINE_FUNCTION
  int histmask(int j) const { return j >> HISTBITS & 1; }
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_NEIGH_HISTORY_KOKKOS_H
#endif // FIX_CLASS
