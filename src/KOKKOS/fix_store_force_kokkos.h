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
FixStyle(store/force/kk,FixStoreForceKokkos<LMPDeviceType>);
FixStyle(store/force/kk/device,FixStoreForceKokkos<LMPDeviceType>);
FixStyle(store/force/kk/host,FixStoreForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_STORE_FORCE_KOKKOS_H
#define LMP_FIX_STORE_FORCE_KOKKOS_H

#include "fix_store_force.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixStoreForce{};

template<class DeviceType>
class FixStoreForceKokkos : public FixStoreForce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixStoreForceKokkos(class LAMMPS *, int, char **);
  ~FixStoreForceKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixStoreForce, const int &) const;

 private:
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_int_1d_randomread d_mask;

  DAT::ttransform_kkfloat_1d_3 k_foriginal;
  typename AT::t_kkfloat_1d_3 d_foriginal;
};

}

#endif
#endif
