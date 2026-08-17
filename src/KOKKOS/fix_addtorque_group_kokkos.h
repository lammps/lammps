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

#ifdef FIX_CLASS
// clang-format off
FixStyle(addtorque/group/kk,FixAddTorqueGroupKokkos<LMPDeviceType>);
FixStyle(addtorque/group/kk/device,FixAddTorqueGroupKokkos<LMPDeviceType>);
FixStyle(addtorque/group/kk/host,FixAddTorqueGroupKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ADDTORQUE_GROUP_KOKKOS_H
#define LMP_FIX_ADDTORQUE_GROUP_KOKKOS_H

#include "fix_addtorque_group.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixAddTorqueGroupMass{};
struct TagFixAddTorqueGroupRmass{};

template<class DeviceType>
class FixAddTorqueGroupKokkos : public FixAddTorqueGroup {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 4;

  FixAddTorqueGroupKokkos(class LAMMPS *, int, char **);
  ~FixAddTorqueGroupKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddTorqueGroupMass, const int &, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddTorqueGroupRmass, const int &, value_type) const;

 private:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  Few<double,3> prd;
  Few<double,6> h;
  int triclinic;

  // set before each kernel launch
  double l_xcm[3];
  double l_omega[3];
  double l_domegadt[3];
  double l_mvv2e;
};

}    // namespace LAMMPS_NS

#endif
#endif
