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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/kk,ComputeTempKokkos<LMPDeviceType>)
ComputeStyle(temp/kk/device,ComputeTempKokkos<LMPDeviceType>)
ComputeStyle(temp/kk/host,ComputeTempKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_TEMP_KOKKOS_H
#define LMP_COMPUTE_TEMP_KOKKOS_H

#include "compute_temp.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

  struct s_CTEMP {
    double t0, t1, t2, t3, t4, t5;
    KOKKOS_INLINE_FUNCTION
    s_CTEMP() {
      t0 = t1 = t2 = t3 = t4 = t5 = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    s_CTEMP& operator+=(const s_CTEMP &rhs){
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

template<int RMASS>
struct TagComputeTempScalar{};

template<int RMASS>
struct TagComputeTempVector{};

template<class DeviceType>
class ComputeTempKokkos : public ComputeTemp {
 public:
  typedef DeviceType device_type;
  typedef CTEMP value_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeTempKokkos(class LAMMPS *, int, char **);
  virtual ~ComputeTempKokkos() {}
  double compute_scalar();
  void compute_vector();

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempScalar<RMASS>, const int&, CTEMP&) const;

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempVector<RMASS>, const int&, CTEMP&) const;

 protected:
  typename ArrayTypes<DeviceType>::t_v_array_randomread v;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread rmass;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread mask;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
