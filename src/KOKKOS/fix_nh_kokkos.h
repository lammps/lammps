// clang-format off
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

#ifndef LMP_FIX_NH_KOKKOS_H
#define LMP_FIX_NH_KOKKOS_H

#include "fix_nh.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int TRICLINIC_FLAG>
struct TagFixNH_nh_v_press{};

template<int RMASS>
struct TagFixNH_nve_v{};

struct TagFixNH_nve_x{};

struct TagFixNH_nh_v_temp{};

template<class DeviceType>
class FixNHKokkos : public FixNH {
 public:
  typedef DeviceType device_type;

  FixNHKokkos(class LAMMPS *, int, char **);

  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void pre_exchange() override;

  template<int TRICLINIC_FLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNH_nh_v_press<TRICLINIC_FLAG>, const int&) const;

  template<int RMASS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNH_nve_v<RMASS>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNH_nve_x, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNH_nh_v_temp, const int&) const;

 protected:
  void remap() override;

  void nve_x() override;            // may be overwritten by child classes
  void nve_v() override;
  void nh_v_press() override;
  void nh_v_temp() override;

  F_FLOAT factor[3];

  class DomainKokkos *domainKK;

  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_f_array_const f;
  typename ArrayTypes<DeviceType>::t_float_1d rmass;
  typename ArrayTypes<DeviceType>::t_float_1d mass;
  typename ArrayTypes<DeviceType>::t_int_1d type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
};

}

#endif

