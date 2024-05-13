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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(third_order/kk,ThirdOrderKokkos);
CommandStyle(third_order/kk/device,ThirdOrderKokkos);
CommandStyle(third_order/kk/host,ThirdOrderKokkos);
// clang-format on
#else

#ifndef LMP_THIRD_ORDER_KOKKOS_H
#define LMP_THIRD_ORDER_KOKKOS_H

#include "kokkos_type.h"
#include "third_order.h"

namespace LAMMPS_NS {

class ThirdOrderKokkos : public ThirdOrder {
 public:
  ThirdOrderKokkos(class LAMMPS *);

  void command(int, char **) override;
  void setup();

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const
  {
    f(i, 0) += f_merge_copy(i, 0);
    f(i, 1) += f_merge_copy(i, 1);
    f(i, 2) += f_merge_copy(i, 2);
  }

 protected:
  void update_force() override;
  void force_clear() override;
  DAT::t_f_array f_merge_copy, f;
};
}    // namespace LAMMPS_NS

#endif    //LMP_THIRD_ORDER_KOKKOS_H
#endif
