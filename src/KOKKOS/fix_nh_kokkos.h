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
  virtual ~FixNHKokkos();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void pre_exchange();

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
  virtual void remap();

  virtual void nve_x();            // may be overwritten by child classes
  virtual void nve_v();
  virtual void nh_v_press();
  virtual void nh_v_temp();

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

/* ERROR/WARNING messages:

E: Cannot (yet) use rigid bodies with fix nh and Kokkos

Self-explanatory.

E: Fix npt/nph has tilted box too far in one step - periodic cell is too far from equilibrium state

Self-explanatory.  The change in the box tilt is too extreme
on a short timescale.

*/
