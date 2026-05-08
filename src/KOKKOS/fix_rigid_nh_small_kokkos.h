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

FixStyle(rigid/nve/small/kk,        FixRigidNVESmallKokkos<LMPDeviceType>);
FixStyle(rigid/nve/small/kk/device, FixRigidNVESmallKokkos<LMPDeviceType>);
FixStyle(rigid/nve/small/kk/host,   FixRigidNVESmallKokkos<LMPHostType>);

FixStyle(rigid/nvt/small/kk,        FixRigidNVTSmallKokkos<LMPDeviceType>);
FixStyle(rigid/nvt/small/kk/device, FixRigidNVTSmallKokkos<LMPDeviceType>);
FixStyle(rigid/nvt/small/kk/host,   FixRigidNVTSmallKokkos<LMPHostType>);

FixStyle(rigid/nph/small/kk,        FixRigidNPHSmallKokkos<LMPDeviceType>);
FixStyle(rigid/nph/small/kk/device, FixRigidNPHSmallKokkos<LMPDeviceType>);
FixStyle(rigid/nph/small/kk/host,   FixRigidNPHSmallKokkos<LMPHostType>);

FixStyle(rigid/npt/small/kk,        FixRigidNPTSmallKokkos<LMPDeviceType>);
FixStyle(rigid/npt/small/kk/device, FixRigidNPTSmallKokkos<LMPDeviceType>);
FixStyle(rigid/npt/small/kk/host,   FixRigidNPTSmallKokkos<LMPHostType>);

// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RIGID_NH_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_NH_SMALL_KOKKOS_H

#include "fix_rigid_nh_small.h"
#include "fix_rigid_base_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType, bool TSTAT, bool PSTAT>
class FixRigidNHSmallKokkos : public FixRigidBaseKokkos<DeviceType, FixRigidNHSmall> {
 public:
  FixRigidNHSmallKokkos(class LAMMPS *, int, char **);

 protected:

  using Pointers::error;

  using Fix::style;

  using FixRigidNHSmall::tstat_flag;
  using FixRigidNHSmall::t_start;
  using FixRigidNHSmall::t_stop;
  using FixRigidNHSmall::t_period;
  using FixRigidNHSmall::pstat_flag;
  using FixRigidNHSmall::p_start;
  using FixRigidNHSmall::p_stop;
  using FixRigidNHSmall::p_flag;
  using FixRigidNHSmall::p_freq;
  using FixRigidNHSmall::p_period;


};

template<class DeviceType>
using FixRigidNVESmallKokkos = FixRigidNHSmallKokkos<DeviceType,false,false>;

template<class DeviceType>
using FixRigidNVTSmallKokkos = FixRigidNHSmallKokkos<DeviceType,true,false>;

template<class DeviceType>
using FixRigidNPHSmallKokkos = FixRigidNHSmallKokkos<DeviceType,false,true>;

template<class DeviceType>
using FixRigidNPTSmallKokkos = FixRigidNHSmallKokkos<DeviceType,true,true>;

}    // namespace LAMMPS_NS

#endif // !LMP_FIX_RIGID_NH_SMALL_KOKKOS_H
#endif
