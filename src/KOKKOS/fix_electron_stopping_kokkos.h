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
FixStyle(electron/stopping/kk,FixElectronStoppingKokkos<LMPDeviceType>);
FixStyle(electron/stopping/kk/device,FixElectronStoppingKokkos<LMPDeviceType>);
FixStyle(electron/stopping/kk/host,FixElectronStoppingKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ELECTRON_STOPPING_KOKKOS_H
#define LMP_FIX_ELECTRON_STOPPING_KOKKOS_H

#include "fix_electron_stopping.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixElectronStopping{};

struct FixElectronStoppingErrorValue {
  int i;
  double energy;
};

KOKKOS_INLINE_FUNCTION bool
operator<(const FixElectronStoppingErrorValue &lhs, const FixElectronStoppingErrorValue &rhs) {
  return lhs.i < rhs.i;
}

template<class DeviceType>
class FixElectronStoppingKokkos : public FixElectronStopping {
 public:

  typedef ArrayTypes<DeviceType> AT;

  FixElectronStoppingKokkos(class LAMMPS *, int, char **);
  ~FixElectronStoppingKokkos() override = default;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixElectronStopping, const int&, double&, FixElectronStoppingErrorValue&) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr_const x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread tag;
  typename AT::t_int_1d_const d_mask;
  typename AT::t_kkfloat_1d_randomread d_mass;
  typename AT::t_kkfloat_1d_const d_rmass;

  typename AT::t_int_1d_const d_numneigh;

  typename AT::t_int_1d d_match;

  typename AT::t_double_2d_const d_elstop_ranges;

  double mvv2e;
  double dt;
};
}

namespace Kokkos { //reduction identity must be defined in Kokkos namespace
   template<>
   struct reduction_identity<LAMMPS_NS::FixElectronStoppingErrorValue> {
      KOKKOS_FORCEINLINE_FUNCTION static LAMMPS_NS::FixElectronStoppingErrorValue min() {
         return {INT_MAX, 0};
      }
   };
}

#endif
#endif
