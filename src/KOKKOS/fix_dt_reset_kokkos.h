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
FixStyle(dt/reset/kk,FixDtResetKokkos<LMPDeviceType>);
FixStyle(dt/reset/kk/device,FixDtResetKokkos<LMPDeviceType>);
FixStyle(dt/reset/kk/host,FixDtResetKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DT_RESET_KOKKOS_H
#define LMP_FIX_DT_RESET_KOKKOS_H

#include "fix_dt_reset.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixDtResetMass{};
struct TagFixDtResetRMass{};

template<class DeviceType>
class FixDtResetKokkos : public FixDtReset {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixDtResetKokkos(class LAMMPS *, int, char **);
//  ~FixDtResetKokkos() override;
  void init() override;
  void end_of_step() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixDtResetMass, const int&, double&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixDtResetRMass, const int&, double&) const;

 private:
  typename AT::t_v_array v;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread rmass;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;


  Kokkos::DualView<double*, Kokkos::LayoutRight, DeviceType> k_emax;
};

}

#endif
#endif

