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
FixStyle(external/kk,FixExternalKokkos<LMPDeviceType>);
FixStyle(external/kk/device,FixExternalKokkos<LMPDeviceType>);
FixStyle(external/kk/host,FixExternalKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EXTERNAL_KOKKOS_H
#define LMP_FIX_EXTERNAL_KOKKOS_H

#include "fix_external.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixExternal{};

template<class DeviceType>
class FixExternalKokkos : public FixExternal {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixExternalKokkos(class LAMMPS *, int, char **);
  ~FixExternalKokkos() override;
  void init() override;
  void post_force(int) override;
  void grow_arrays(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixExternal, const int&) const;

 private:
  DAT::ttransform_kkfloat_2d k_fexternal;
  typename AT::t_kkfloat_2d_randomread d_fexternal;

  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkacc_1d_3 f;
};

}

#endif
#endif
