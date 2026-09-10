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
ComputeStyle(temp/region/kk,ComputeTempRegionKokkos<LMPDeviceType>);
ComputeStyle(temp/region/kk/device,ComputeTempRegionKokkos<LMPDeviceType>);
ComputeStyle(temp/region/kk/host,ComputeTempRegionKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_REGION_KOKKOS_H
#define LMP_COMPUTE_TEMP_REGION_KOKKOS_H

#include "compute_temp_region.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempRegionScalar{};

template<int RMASS>
struct TagComputeTempRegionVector{};

struct TagComputeTempRegionRemoveBias{};

struct TagComputeTempRegionRestoreBias{};

template<class DeviceType>
class ComputeTempRegionKokkos : public ComputeTempRegion {
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

  ComputeTempRegionKokkos(class LAMMPS *, int, char **);
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;

  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;
  void restore_bias_all_kk() override;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRegionScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRegionVector<RMASS>, const int&, CTEMP&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRegionRemoveBias, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempRegionRestoreBias, const int&) const;

 protected:
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 vbiasall;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;

  DAT::tdual_int_1d k_match;
  typename AT::t_int_1d_randomread d_match;

  // evaluate the region for all atoms in the group into d_match

  void region_match_all();
};

}    // namespace LAMMPS_NS

#endif
#endif
