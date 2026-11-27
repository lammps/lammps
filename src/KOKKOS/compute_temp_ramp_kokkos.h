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
ComputeStyle(temp/ramp/kk,ComputeTempRampKokkos<LMPDeviceType>);
ComputeStyle(temp/ramp/kk/device,ComputeTempRampKokkos<LMPDeviceType>);
ComputeStyle(temp/ramp/kk/host,ComputeTempRampKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_RAMP_KOKKOS_H
#define LMP_COMPUTE_TEMP_RAMP_KOKKOS_H

#include "compute_temp_ramp.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempRampScalar{};

template<int RMASS>
struct TagComputeTempRampVector{};

struct TagComputeTempRampRemoveBias{};

struct TagComputeTempRampRestoreBias{};

template<class DeviceType>
class ComputeTempRampKokkos : public ComputeTempRamp {
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
      t0 += rhs.t0; t1 += rhs.t1; t2 += rhs.t2;
      t3 += rhs.t3; t4 += rhs.t4; t5 += rhs.t5;
      return *this;
    }
  };

  typedef s_CTEMP CTEMP;
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef CTEMP value_type;

  ComputeTempRampKokkos(class LAMMPS *, int, char **);
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;
  void restore_bias_all_kk() override;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRampScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRampVector<RMASS>, const int&, CTEMP&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRampRemoveBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRampRestoreBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  double ramp_bias(const int &i) const {
    double fraction = (static_cast<double>(x(i,coord_dim)) - coord_lo) / (coord_hi - coord_lo);
    fraction = (fraction < 0.0) ? 0.0 : fraction;
    fraction = (fraction > 1.0) ? 1.0 : fraction;
    return v_lo + fraction*(v_hi - v_lo);
  }

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
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
