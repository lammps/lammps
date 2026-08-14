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
FixStyle(nve/noforce/kk,FixNVENoforceKokkos<LMPDeviceType>);
FixStyle(nve/noforce/kk/device,FixNVENoforceKokkos<LMPDeviceType>);
FixStyle(nve/noforce/kk/host,FixNVENoforceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NVE_NOFORCE_KOKKOS_H
#define LMP_FIX_NVE_NOFORCE_KOKKOS_H

#include "fix_nve_noforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixNVENoforce{};

template<class DeviceType>
class FixNVENoforceKokkos : public FixNVENoforce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNVENoforceKokkos(class LAMMPS *, int, char **);
  ~FixNVENoforceKokkos() override;
  void init() override;
  void initial_integrate(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVENoforce, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_int_1d_randomread mask;

};

}

#endif
#endif
