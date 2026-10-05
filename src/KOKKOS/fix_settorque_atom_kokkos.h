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
FixStyle(settorque/atom/kk,FixSetTorqueAtomKokkos<LMPDeviceType>);
FixStyle(settorque/atom/kk/device,FixSetTorqueAtomKokkos<LMPDeviceType>);
FixStyle(settorque/atom/kk/host,FixSetTorqueAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SETTORQUE_ATOM_KOKKOS_H
#define LMP_FIX_SETTORQUE_ATOM_KOKKOS_H

#include "fix_settorque_atom.h"

#include "fix_setforce_kokkos.h"    // for double_3
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixSetTorqueAtomConstant{};
struct TagFixSetTorqueAtomNonConstant{};

template<class DeviceType>
class FixSetTorqueAtomKokkos : public FixSetTorqueAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double_3 value_type;

  FixSetTorqueAtomKokkos(class LAMMPS *, int, char **);
  ~FixSetTorqueAtomKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetTorqueAtomConstant, const int&, double_3&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSetTorqueAtomNonConstant, const int&, double_3&) const;

 private:
  DAT::ttransform_kkfloat_2d k_storque;
  typename AT::t_kkfloat_2d_randomread d_storque;
  typename AT::t_int_1d d_match;

  typename AT::t_kkacc_1d_3 d_torque;
  typename AT::t_int_1d_randomread d_mask;
};

}

#endif
#endif
