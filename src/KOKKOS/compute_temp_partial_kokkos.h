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
ComputeStyle(temp/partial/kk,ComputeTempPartialKokkos<LMPDeviceType>);
ComputeStyle(temp/partial/kk/device,ComputeTempPartialKokkos<LMPDeviceType>);
ComputeStyle(temp/partial/kk/host,ComputeTempPartialKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_PARTIAL_KOKKOS_H
#define LMP_COMPUTE_TEMP_PARTIAL_KOKKOS_H

#include "compute_temp_partial.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempPartialScalar{};

template<int RMASS>
struct TagComputeTempPartialVector{};

struct TagComputeTempPartialRemoveBias{};

struct TagComputeTempPartialRestoreBias{};

struct TagComputeTempPartialReapplyBias{};

template<class DeviceType>
class ComputeTempPartialKokkos : public ComputeTempPartial {
 public:

  struct s_CTEMP {
    double t0, t1, t2, t3, t4, t5;
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    s_CTEMP() {
      t0 = t1 = t2 = t3 = t4 = t5 = 0.0;
    }
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    s_CTEMP& operator+=(const s_CTEMP &rhs) {
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

  ComputeTempPartialKokkos(class LAMMPS *, int, char **);
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;
  void restore_bias_all_kk() override;
  void reapply_bias_all() override;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempPartialScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempPartialVector<RMASS>, const int&, CTEMP&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempPartialRemoveBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempPartialRestoreBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempPartialReapplyBias, const int&) const;

 protected:
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 vbiasall;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
};

}    // namespace LAMMPS_NS

#endif
#endif
