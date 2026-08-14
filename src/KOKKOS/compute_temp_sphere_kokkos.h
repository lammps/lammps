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

// Tags for Kokkos kernels

template <int MODE>
struct TagComputeTempSphereScalar {};

template <int MODE>
struct TagComputeTempSphereVector {};

template <class DeviceType>
class ComputeTempSphereKokkos : public ComputeTempSphere {
 public:
  struct s_CTEMP {
    double t0, t1, t2, t3, t4, t5;
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    s_CTEMP() { t0 = t1 = t2 = t3 = t4 = t5 = 0.0; }
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    s_CTEMP &operator+=(const s_CTEMP &rhs)
    {
      t0 += rhs.t0;
      t1 += rhs.t1;
      t2 += rhs.t2;
      t3 += rhs.t3;
      t4 += rhs.t4;
      t5 += rhs.t5;
      return *this;
    }
  };

  typedef s_CTEMP CTEMP;
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef CTEMP value_type;

  ComputeTempSphereKokkos(class LAMMPS *, int, char **);

  double compute_scalar() override;
  void compute_vector() override;

  template <int MODE>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempSphereScalar<MODE>, const int &, CTEMP &) const;

  template <int MODE>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempSphereVector<MODE>, const int &, CTEMP &) const;

 private:
  class AtomKokkos *atomKK;

  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkfloat_1d_3_randomread omega_kk;
  typename AT::t_kkfloat_1d_randomread rmass_kk;
  typename AT::t_kkfloat_1d_randomread radius_kk;
  typename AT::t_int_1d_randomread mask;
};

}    // namespace LAMMPS_NS

#endif
#endif
