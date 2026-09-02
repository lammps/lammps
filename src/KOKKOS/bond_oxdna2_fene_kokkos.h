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

#ifdef BOND_CLASS
// clang-format off
BondStyle(oxdna2/fene/kk,BondOxdna2FENEKokkos<LMPDeviceType>);
BondStyle(oxdna2/fene/kk/device,BondOxdna2FENEKokkos<LMPDeviceType>);
BondStyle(oxdna2/fene/kk/host,BondOxdna2FENEKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_BOND_OXDNA2_FENE_KOKKOS_H
#define LMP_BOND_OXDNA2_FENE_KOKKOS_H

#include "bond_oxdna_fene_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class BondOxdna2FENEKokkos : public BondOxdnaFENEKokkos<DeviceType> {
 public:
  BondOxdna2FENEKokkos(class LAMMPS *);
  ~BondOxdna2FENEKokkos() {}
};

template<class DeviceType>
BondOxdna2FENEKokkos<DeviceType>::BondOxdna2FENEKokkos(LAMMPS *lmp) : BondOxdnaFENEKokkos<DeviceType>(lmp)
{
   this->oxdnaflag = BondOxdnaFENEKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA2;
}
}    // namespace LAMMPS_NS

#endif
#endif
