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
ComputeStyle(temp/profile/kk,ComputeTempProfileKokkos<LMPDeviceType>);
ComputeStyle(temp/profile/kk/device,ComputeTempProfileKokkos<LMPDeviceType>);
ComputeStyle(temp/profile/kk/host,ComputeTempProfileKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_PROFILE_KOKKOS_H
#define LMP_COMPUTE_TEMP_PROFILE_KOKKOS_H

#include "compute_temp_profile.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeTempProfileBin{};

template<int RMASS>
struct TagComputeTempProfileScatter{};

template<int RMASS>
struct TagComputeTempProfileScalar{};

template<int RMASS>
struct TagComputeTempProfileVector{};

struct TagComputeTempProfileRemoveBias{};

struct TagComputeTempProfileRestoreBias{};

template<class DeviceType>
class ComputeTempProfileKokkos : public ComputeTempProfile {
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

  ComputeTempProfileKokkos(class LAMMPS *, int, char **);

  double compute_scalar() override;
  void compute_vector() override;
  void compute_array() override;

  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;
  void restore_bias_all_kk() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileBin, const int&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileScatter<RMASS>, const int&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileVector<RMASS>, const int&, CTEMP&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileRemoveBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempProfileRestoreBias, const int&) const;

 protected:
  // average COM velocity per bin (device), mirroring ComputeTempProfile::bin_average
  void bin_average_kk();

  typedef Kokkos::View<double**, DeviceType> t_double_2d;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;

  t_double_2d d_vbin;                    // per-bin mass-weighted v + mass + count (scatter)
  t_double_2d d_binave;                  // per-bin COM velocity (after Allreduce + divide)
  typename t_double_2d::host_mirror_type h_vbin, h_binave;   // persistent host mirrors
  typename AT::t_int_1d d_bin;           // per-atom bin index

  int maxbin;

  // binning frame copied to device-friendly scalars (orthogonal path)
  KK_FLOAT m_boxlo[3], m_boxhi[3], m_prd[3], m_invdelta[3];
  int m_periodicity[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
