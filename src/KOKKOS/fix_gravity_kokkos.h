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

FixStyle(gravity/kk,FixGravityKokkos<Device>)
FixStyle(gravity/kk/device,FixGravityKokkos<Device>)
FixStyle(gravity/kk/host,FixGravityKokkos<Host>)

#else

#ifndef LMP_FIX_GRAVITY_KOKKOS_H
#define LMP_FIX_GRAVITY_KOKKOS_H

#include "fix_gravity.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixGravityRMass {};
struct TagFixGravityMass {};

template<ExecutionSpace Space>
class FixGravityKokkos : public FixGravity {
  public:
    typedef typename GetDeviceType<Space>::value DeviceType;
    typedef ArrayTypes<Space> AT;

    FixGravityKokkos(class LAMMPS *, int, char **);
    virtual ~FixGravityKokkos() {}
    void post_force(int);

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixGravityRMass, const int, KK_FLOAT &) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixGravityMass, const int, KK_FLOAT &) const;

  private:
    typename AT::t_float_1d_3_lr x;
    typename AT::t_float_1d_3 f;
    typename AT::t_float_1d_randomread rmass;
    typename AT::t_float_1d_randomread mass;
    typename AT::t_int_1d type;
    typename AT::t_int_1d mask;
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_GRAVITY_KOKKOS_H
#endif // FIX_CLASS
