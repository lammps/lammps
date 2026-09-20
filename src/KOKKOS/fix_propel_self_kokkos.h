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
FixStyle(propel/self/kk,FixPropelSelfKokkos<LMPDeviceType>);
FixStyle(propel/self/kk/device,FixPropelSelfKokkos<LMPDeviceType>);
FixStyle(propel/self/kk/host,FixPropelSelfKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_PROPEL_SELF_KOKKOS_H
#define LMP_FIX_PROPEL_SELF_KOKKOS_H

#include "fix_propel_self.h"

#include "kokkos_few.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixPropelSelfDipole{};
struct TagFixPropelSelfVelocity{};

template<class DeviceType>
class FixPropelSelfKokkos : public FixPropelSelf {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 6;

  FixPropelSelfKokkos(class LAMMPS *, int, char **);
  ~FixPropelSelfKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixPropelSelfDipole, const int&, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixPropelSelfVelocity, const int&, value_type) const;

 private:
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_kkfloat_1d_3_randomread d_v;
  typename AT::t_kkfloat_1d_4_randomread d_mu;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_imageint_1d_randomread d_image;
  typename AT::t_int_1d_randomread d_mask;

  Few<double,3> prd;
  Few<double,6> h;
  int triclinic;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void tally(value_type, const int &, KK_FLOAT, KK_FLOAT, KK_FLOAT) const;
};

}

#endif
#endif
