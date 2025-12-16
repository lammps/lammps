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
ComputeStyle(temp/com/kk,ComputeTempCOMKokkos<LMPDeviceType>);
ComputeStyle(temp/com/kk/device,ComputeTempCOMKokkos<LMPDeviceType>);
ComputeStyle(temp/com/kk/host,ComputeTempCOMKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_COM_KOKKOS_H
#define LMP_COMPUTE_TEMP_COM_KOKKOS_H

#include "compute_temp_com.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempCOMScalar{};

template<int RMASS>
struct TagComputeTempCOMVector{};

struct TagComputeTempCOMRemoveBias{};
struct TagComputeTempCOMRestoreBias{};

template<class DeviceType>
class ComputeTempCOMKokkos : public ComputeTempCOM {
 public:

  struct s_CTEMP {
    double t0, t1, t2, t3, t4, t5;
    KOKKOS_INLINE_FUNCTION
    s_CTEMP() {
      t0 = t1 = t2 = t3 = t4 = t5 = 0.0;
    }
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

  ComputeTempCOMKokkos(class LAMMPS *, int, char **);
  double compute_scalar() override;
  void compute_vector() override;
  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempCOMScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempCOMVector<RMASS>, const int&, CTEMP&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempCOMRemoveBias, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempCOMRestoreBias, const int &i) const;

 protected:
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;

  class GroupKokkos *groupKK;
};

}    // namespace LAMMPS_NS

#endif
#endif
