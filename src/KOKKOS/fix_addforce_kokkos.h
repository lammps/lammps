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
FixStyle(addforce/kk,FixAddForceKokkos<LMPDeviceType>);
FixStyle(addforce/kk/device,FixAddForceKokkos<LMPDeviceType>);
FixStyle(addforce/kk/host,FixAddForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ADD_FORCE_KOKKOS_H
#define LMP_FIX_ADD_FORCE_KOKKOS_H

#include "fix_addforce.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

struct TagFixAddForceConstant{};
struct TagFixAddForceNonConstant{};

template<class DeviceType>
class FixAddForceKokkos : public FixAddForce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 10;

  FixAddForceKokkos(class LAMMPS *, int, char **);
  ~FixAddForceKokkos() override;
  void init() override;
  void post_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceConstant, const int&, value_type) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceNonConstant, const int&, value_type) const;

 private:
  DAT::ttransform_kkfloat_2d k_sforce;
  typename AT::t_kkfloat_2d_randomread d_sforce;
  typename AT::t_int_1d d_match;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;

  Few<double,3> prd;
  Few<double,6> h;
  int triclinic;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  KOKKOS_INLINE_FUNCTION
  void v_tally(value_type, int, KK_FLOAT*) const;
};

}

#endif
#endif
