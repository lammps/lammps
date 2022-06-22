/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS
// clang-format off
RegionStyle(block/kk,RegBlockKokkos<LMPDeviceType>);
RegionStyle(block/kk/device,RegBlockKokkos<LMPDeviceType>);
RegionStyle(block/kk/host,RegBlockKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_BLOCK_KOKKOS_H
#define LMP_REGION_BLOCK_KOKKOS_H

#include "region_block.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagRegBlockMatchAll{};

template<class DeviceType>
class RegBlockKokkos : public RegBlock, public KokkosBase {
  friend class FixPour;

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegBlockKokkos(class LAMMPS *, int, char **);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegBlockMatchAll, const int&) const;

 private:
  int groupbit;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_int_1d_randomread mask;

  KOKKOS_INLINE_FUNCTION
  int k_inside(double, double, double) const;
  KOKKOS_INLINE_FUNCTION
  int match(double, double, double) const;
  KOKKOS_INLINE_FUNCTION
  void inverse_transform(double &, double &, double &) const;
  KOKKOS_INLINE_FUNCTION
  void rotate(double &, double &, double &, double) const;

};

}

#endif
#endif

