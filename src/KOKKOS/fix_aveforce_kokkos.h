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
FixStyle(aveforce/kk,FixAveForceKokkos<LMPDeviceType>);
FixStyle(aveforce/kk/device,FixAveForceKokkos<LMPDeviceType>);
FixStyle(aveforce/kk/host,FixAveForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_AVEFORCE_KOKKOS_H
#define LMP_FIX_AVEFORCE_KOKKOS_H

#include "fix_aveforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixAveForceReduce{};
struct TagFixAveForceApply{};

template<class DeviceType>
class FixAveForceKokkos : public FixAveForce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 4;

  FixAveForceKokkos(class LAMMPS *, int, char **);
  ~FixAveForceKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveForceReduce, const int &, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveForceApply, const int &) const;

 private:
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d d_match;

  double m_fave[3];

};

}

#endif
#endif
