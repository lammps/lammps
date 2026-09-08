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
FixStyle(oneway/kk,FixOneWayKokkos<LMPDeviceType>);
FixStyle(oneway/kk/device,FixOneWayKokkos<LMPDeviceType>);
FixStyle(oneway/kk/host,FixOneWayKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ONEWAY_KOKKOS_H
#define LMP_FIX_ONEWAY_KOKKOS_H

#include "fix_oneway.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixOneWay{};

template<class DeviceType>
class FixOneWayKokkos : public FixOneWay {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixOneWayKokkos(class LAMMPS *, int, char **);
  ~FixOneWayKokkos() override;
  void init() override;
  void end_of_step() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOneWay, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d d_match;

};

}

#endif
#endif
