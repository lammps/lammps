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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/sphere/kk,ComputeTempSphereKokkos<LMPDeviceType>);
ComputeStyle(temp/sphere/kk/device,ComputeTempSphereKokkos<LMPDeviceType>);
ComputeStyle(temp/sphere/kk/host,ComputeTempSphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_SPHERE_KOKKOS_H
#define LMP_COMPUTE_TEMP_SPHERE_KOKKOS_H

#include "compute_temp_sphere.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeTempSphereScalar{};
struct TagComputeTempSphereVector{};

template<class DeviceType>
class ComputeTempSphereKokkos : public ComputeTempSphere {
 public:

  struct s_CTEMP {
    double t0, t1, t2, t3, t4, t5;
    KOKKOS_INLINE_FUNCTION
    s_CTEMP() { t0 = t1 = t2 = t3 = t4 = t5 = 0.0; }
    KOKKOS_INLINE_FUNCTION
    s_CTEMP& operator+=(const s_CTEMP &rhs) {
      t0 += rhs.t0; t1 += rhs.t1; t2 += rhs.t2;
      t3 += rhs.t3; t4 += rhs.t4; t5 += rhs.t5;
      return *this;
    }
  };

  typedef s_CTEMP CTEMP;
  typedef DeviceType device_type;
  typedef CTEMP value_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeTempSphereKokkos(class LAMMPS *, int, char **);

  double compute_scalar() override;
  void compute_vector() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempSphereScalar, const int &, CTEMP &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempSphereVector, const int &, CTEMP &) const;

 protected:
  class AtomKokkos *atomKK;
  typename AT::t_v_array_randomread v;
  typename AT::t_v_array_randomread omega;
  typename AT::t_float_1d_randomread radius;
  typename AT::t_float_1d_randomread rmass;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif
