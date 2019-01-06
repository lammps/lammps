/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gravity/kk,FixGravityKokkos<LMPDeviceType>)
FixStyle(gravity/kk/device,FixGravityKokkos<LMPDeviceType>)
FixStyle(gravity/kk/host,FixGravityKokkos<LMPHostType>)

#else

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
    virtual ~FixGravityKokkos() {}
    void post_force(int);

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
