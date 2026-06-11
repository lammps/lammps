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
FixStyle(halt/kk,FixHaltKokkos<LMPDeviceType>);
FixStyle(halt/kk/device,FixHaltKokkos<LMPDeviceType>);
FixStyle(halt/kk/host,FixHaltKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_HALT_KOKKOS_H
#define LMP_FIX_HALT_KOKKOS_H

#include "fix_halt.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixHaltBondmax{};

template<class DeviceType>
class FixHaltKokkos : public FixHalt {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixHaltKokkos(class LAMMPS *, int, char **);
  ~FixHaltKokkos() {};
  void end_of_step() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixHaltBondmax, const int &, KK_FLOAT &) const;

 protected:
  double bondmax() override;

 private:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_2d_lr d_bondlist;

  class NeighborKokkos *neighborKK;
};

}

#endif
#endif
