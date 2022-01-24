/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/deform/kk,ComputeTempDeformKokkos<LMPDeviceType>);
ComputeStyle(temp/deform/kk/device,ComputeTempDeformKokkos<LMPDeviceType>);
ComputeStyle(temp/deform/kk/host,ComputeTempDeformKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_DEFORM_KOKKOS_H
#define LMP_COMPUTE_TEMP_DEFORM_KOKKOS_H

#include "compute_temp_deform.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempDeformScalar{};

template<int RMASS>
struct TagComputeTempDeformVector{};

struct TagComputeTempDeformRemoveBias{};

struct TagComputeTempDeformRestoreBias{};

template<class DeviceType>
class ComputeTempDeformKokkos: public ComputeTempDeform {
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

    KOKKOS_INLINE_FUNCTION
    void operator+=(const volatile s_CTEMP &rhs) volatile {
      t0 += rhs.t0;
      t1 += rhs.t1;
      t2 += rhs.t2;
      t3 += rhs.t3;
      t4 += rhs.t4;
      t5 += rhs.t5;
    }
  };

  typedef s_CTEMP CTEMP;
  typedef CTEMP value_type;
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeTempDeformKokkos(class LAMMPS *, int, char **);
  ~ComputeTempDeformKokkos();
  double compute_scalar();
  void compute_vector();
  void remove_bias_all();
  void restore_bias_all();

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformVector<RMASS>, const int&, CTEMP&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformRemoveBias, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformRestoreBias, const int &i) const;

 protected:
  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_v_array vbiasall;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread rmass;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread mask;

  class DomainKokkos *domainKK;

  Few<double, 6> h_rate, h_ratelo;

  };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
