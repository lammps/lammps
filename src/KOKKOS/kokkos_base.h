// clang-format off
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

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"

#include <Kokkos_Sort.hpp>

namespace LAMMPS_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  // Pair
  virtual int pack_forward_comm_kokkos(int, DAT::tdual_int_2d,
                                       int, DAT::tdual_xfloat_1d &,
                                       int, int *) {return 0;};
  virtual void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) {}

  virtual int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) {return 0;};
  virtual void unpack_reverse_comm_kokkos(int, DAT::tdual_int_2d,
                                          int, DAT::tdual_xfloat_1d &) {}

  // Fix
  virtual int pack_forward_comm_fix_kokkos(int, DAT::tdual_int_2d,
                                           int, DAT::tdual_xfloat_1d &,
                                           int, int *) {return 0;};
  virtual void unpack_forward_comm_fix_kokkos(int, int, DAT::tdual_xfloat_1d &) {}

  virtual int pack_exchange_kokkos(const int & /*nsend*/, DAT::tdual_xfloat_2d & /*k_buf*/,
                                   DAT::tdual_int_1d /*k_sendlist*/,
                                   DAT::tdual_int_1d /*k_copylist*/,
                                   ExecutionSpace /*space*/) { return 0; }
  virtual void unpack_exchange_kokkos(DAT::tdual_xfloat_2d & /*k_buf*/,
                                      DAT::tdual_int_1d & /*indices*/, int /*nrecv*/,
                                      ExecutionSpace /*space*/) {}

  // Region
  virtual void match_all_kokkos(int, DAT::tdual_int_1d) {}

  using KeyViewType = DAT::t_x_array;
  using BinOp = BinOp3DLAMMPS<KeyViewType>;
  virtual void
    sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> & /*Sorter*/) {}
};

}

#endif

