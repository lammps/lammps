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
FixStyle(viscous/nonlinear/kk,FixViscousNonlinearKokkos<LMPDeviceType>);
FixStyle(viscous/nonlinear/kk/device,FixViscousNonlinearKokkos<LMPDeviceType>);
FixStyle(viscous/nonlinear/kk/host,FixViscousNonlinearKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_VISCOUS_NONLINEAR_KOKKOS_H
#define LMP_FIX_VISCOUS_NONLINEAR_KOKKOS_H

#include "fix_viscous_nonlinear.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixViscousNonlinear{};

template<class DeviceType>
class FixViscousNonlinearKokkos : public FixViscousNonlinear {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixViscousNonlinearKokkos(class LAMMPS *, int, char **);
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixViscousNonlinear, const int &) const;

 private:
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d_3_randomread d_v;
  typename AT::t_kkfloat_1d_randomread d_radius;
  typename AT::t_int_1d_randomread d_mask;

  KK_FLOAT m_v_fluid[3];
};

}

#endif
#endif
