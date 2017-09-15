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

FixStyle(setforce/kk,FixSetForceKokkos<LMPDeviceType>)
FixStyle(setforce/kk/device,FixSetForceKokkos<LMPDeviceType>)
FixStyle(setforce/kk/host,FixSetForceKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_SET_FORCE_KOKKOS_H
#define LMP_FIX_SET_FORCE_KOKKOS_H

#include "fix_setforce.h"
#include "region.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct s_double_3 {
  double d0, d1, d2;
  KOKKOS_INLINE_FUNCTION
  s_double_3() {
    d0 = d1 = d2 = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  s_double_3& operator+=(const s_double_3 &rhs){
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_double_3 &rhs) volatile {
    d0 += rhs.d0;
    d1 += rhs.d1;
    d2 += rhs.d2;
  }
};
typedef s_double_3 double_3;

struct TagFixSetForceConstant{};

struct TagFixSetForceNonConstant{};

template<class DeviceType>
class FixSetForceKokkos : public FixSetForce {
 public:
  typedef DeviceType device_type;
  typedef double_3 value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixSetForceKokkos(class LAMMPS *, int, char **);
  ~FixSetForceKokkos();
  void init();
  void post_force(int);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetForceConstant, const int&, double_3&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetForceNonConstant, const int&, double_3&) const;

 private:
  DAT::tdual_ffloat_2d k_sforce;
  typename AT::t_ffloat_2d_randomread d_sforce;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;

  Region* region;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot (yet) use respa with Kokkos

Self-explanatory.

*/
