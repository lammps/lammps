/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(NEIGH_HISTORY/KK,FixNeighHistoryKokkos<Device>)
FixStyle(NEIGH_HISTORY/KK/DEVICE,FixNeighHistoryKokkos<Device>)
FixStyle(NEIGH_HISTORY/KK/HOST,FixNeighHistoryKokkos<Host>)

#else

#ifndef LMP_FIX_NEIGH_HISTORY_KOKKOS_H
#define LMP_FIX_NEIGH_HISTORY_KOKKOS_H

#include "fix_neigh_history.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {
template <ExecutionSpace Space>
class FixNeighHistoryKokkos : public FixNeighHistory {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef ArrayTypes<Space> AT;

  FixNeighHistoryKokkos(class LAMMPS *, int, char **);
  ~FixNeighHistoryKokkos();

  void init();
  void pre_exchange();
  void setup_post_neighbor();
  virtual void post_neighbor();
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  KOKKOS_INLINE_FUNCTION
  void zero_partner_count_item(const int &i) const;
  KOKKOS_INLINE_FUNCTION
  void pre_exchange_item(const int &ii) const;
  KOKKOS_INLINE_FUNCTION
  void post_neighbor_item(const int &ii) const;

  typename Kokkos::View<int**> d_firstflag;
  typename Kokkos::View<typename AT::t_float_1d::data_type*> d_firstvalue;

 private:
  DAT::tdual_int_1d k_npartner;
  DAT::tdual_tagint_2d k_partner;
  DAT::tdual_float_2d k_valuepartner;

  // for neighbor list lookup
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d d_npartner;
  typename AT::t_tagint_2d d_partner;
  typename AT::t_float_2d d_valuepartner;

  typename AT::t_int_scalar d_resize;
  HAT::t_int_scalar h_resize;
};

template <ExecutionSpace Space>
struct FixNeighHistoryKokkosZeroPartnerCountFunctor {
  FixNeighHistoryKokkos<Space> c;
  FixNeighHistoryKokkosZeroPartnerCountFunctor(FixNeighHistoryKokkos<Space> *c_ptr): c(*c_ptr) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    c.zero_partner_count_item(i);
  }
};

template <ExecutionSpace Space>
struct FixNeighHistoryKokkosPreExchangeFunctor {
  FixNeighHistoryKokkos<Space> c;
  FixNeighHistoryKokkosPreExchangeFunctor(FixNeighHistoryKokkos<Space> *c_ptr): c(*c_ptr) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    c.pre_exchange_item(i);
  }
};

template <ExecutionSpace Space>
struct FixNeighHistoryKokkosPostNeighborFunctor {
  FixNeighHistoryKokkos<Space> c;
  FixNeighHistoryKokkosPostNeighborFunctor(FixNeighHistoryKokkos<Space> *c_ptr): c(*c_ptr) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    c.post_neighbor_item(i);
  }
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_NEIGH_HISTORY_KOKKOS_H
#endif // FIX_CLASS
