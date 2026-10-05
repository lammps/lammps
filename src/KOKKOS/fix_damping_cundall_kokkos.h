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
FixStyle(damping/cundall/kk,FixDampingCundallKokkos<LMPDeviceType>);
FixStyle(damping/cundall/kk/device,FixDampingCundallKokkos<LMPDeviceType>);
FixStyle(damping/cundall/kk/host,FixDampingCundallKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DAMPING_CUNDALL_KOKKOS_H
#define LMP_FIX_DAMPING_CUNDALL_KOKKOS_H

#include "fix_damping_cundall.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

// SCALE selects where the per-atom damping prefactor comes from and matches
// the NONE/TYPE/VARIABLE enum of the base style

template<int SCALE>
struct TagFixDampingCundall{};

template<class DeviceType>
class FixDampingCundallKokkos : public FixDampingCundall {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixDampingCundallKokkos(class LAMMPS *, int, char **);
  ~FixDampingCundallKokkos() override;
  void init() override;
  void post_force(int) override;

  template<int SCALE>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixDampingCundall<SCALE>, const int &) const;

 private:
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkacc_1d_3 d_torque;
  typename AT::t_kkfloat_1d_3_randomread d_v;
  typename AT::t_kkfloat_1d_3_randomread d_omega;
  typename AT::t_int_1d_randomread d_mask;
  typename AT::t_int_1d_randomread d_type;

  // per-type scale factors, and the host-evaluated atom-style variable

  DAT::tdual_kkfloat_1d k_scalegamma;
  typename AT::t_kkfloat_1d_randomread d_scalegamma;
  DAT::tdual_kkfloat_1d k_scaleval;
  typename AT::t_kkfloat_1d_randomread d_scaleval;
  int maxatom_scaleval;
};

}

#endif
#endif
