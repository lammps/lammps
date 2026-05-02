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

FixStyle(efield/kk,        FixEfieldKokkos<LMPDeviceType,false>);
FixStyle(efield/kk/device, FixEfieldKokkos<LMPDeviceType,false>);
FixStyle(efield/kk/host,   FixEfieldKokkos<LMPHostType,false>);

FixStyle(efield/tip4p/kk,        FixEfieldKokkos<LMPDeviceType,true>);
FixStyle(efield/tip4p/kk/device, FixEfieldKokkos<LMPDeviceType,true>);
FixStyle(efield/tip4p/kk/host,   FixEfieldKokkos<LMPHostType,true>);

// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EFIELD_KOKKOS_H
#define LMP_FIX_EFIELD_KOKKOS_H

#include "fix_efield.h"
#include "fix_efield_tip4p.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

template<class DeviceType, bool TIP4P>
class FixEfieldKokkos :
#ifdef LMP_EXTRA_FIX
public std::conditional_t<TIP4P, FixEfieldTIP4P, FixEfield>
#else
public FixEfield
#endif
{
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixEfieldKokkos(class LAMMPS *, int, char **);
  ~FixEfieldKokkos() override;
  void init() override;
  void post_force(int) override;

 private:

  DAT::ttransform_kkfloat_2d k_efield;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;


};


//template <class DeviceType>
//using FixEfieldPlainKokkos = FixEfieldKokkos<DeviceType, false>;

//template <class DeviceType>
//using FixEfieldTIP4PKokkos = FixEfieldKokkos<DeviceType, true>;

} // namespace LAMMPS_NS

#endif // !LMP_FIX_EFIELD_KOKKOS_H
#endif

