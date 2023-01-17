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
FixStyle(gravity/kk,FixGravityKokkos<LMPDeviceType>);
FixStyle(gravity/kk/device,FixGravityKokkos<LMPDeviceType>);
FixStyle(gravity/kk/host,FixGravityKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_GRAVITY_KOKKOS_H
#define LMP_FIX_GRAVITY_KOKKOS_H

#include "fix_gravity.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixGravityRMass {};
struct TagFixGravityMass {};

template<class DeviceType>
class FixGravityKokkos : public FixGravity {
  public:
    FixGravityKokkos(class LAMMPS *, int, char **);

    void post_force(int) override;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixGravityRMass, const int, double &) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixGravityMass, const int, double &) const;

  private:
    typename ArrayTypes<DeviceType>::t_x_array x;
    typename ArrayTypes<DeviceType>::t_f_array f;
    typename ArrayTypes<DeviceType>::t_float_1d_randomread rmass;
    typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
    typename ArrayTypes<DeviceType>::t_int_1d type;
    typename ArrayTypes<DeviceType>::t_int_1d mask;
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_GRAVITY_KOKKOS_H
#endif // FIX_CLASS
