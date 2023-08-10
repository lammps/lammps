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
FixStyle(spring/self/kk,FixSpringSelfKokkos<LMPDeviceType>);
FixStyle(spring/self/kk/device,FixSpringSelfKokkos<LMPDeviceType>);
FixStyle(spring/self/kk/host,FixSpringSelfKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SPRING_SELF_KOKKOS_H
#define LMP_FIX_SPRING_SELF_KOKKOS_H

#include "fix_spring_self.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixSpringSelfKokkos : public FixSpringSelf {
 public:
  typedef DeviceType device_type;
  typedef double value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixSpringSelfKokkos(class LAMMPS *, int, char **);
  ~FixSpringSelfKokkos() override;
  void init() override;
  void post_force(int) override;

 private:
  DAT::tdual_ffloat_2d k_xoriginal;
  typename AT::t_ffloat_2d_randomread d_xoriginal;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif

