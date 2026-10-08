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
BondStyle(oxrna2/fene/kk,BondOxrna2FENEKokkos<LMPDeviceType>);
BondStyle(oxrna2/fene/kk/device,BondOxrna2FENEKokkos<LMPDeviceType>);
BondStyle(oxrna2/fene/kk/host,BondOxrna2FENEKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_BOND_OXRNA2_FENE_KOKKOS_H
#define LMP_BOND_OXRNA2_FENE_KOKKOS_H

#include "bond_oxdna_fene_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class BondOxrna2FENEKokkos : public BondOxdnaFENEKokkos<DeviceType> {
 public:
  BondOxrna2FENEKokkos(class LAMMPS *);
  ~BondOxrna2FENEKokkos() {}
};

template<class DeviceType>
BondOxrna2FENEKokkos<DeviceType>::BondOxrna2FENEKokkos(LAMMPS *lmp) : BondOxdnaFENEKokkos<DeviceType>(lmp)
{
   this->oxdnaflag = BondOxdnaFENEKokkos<DeviceType>::EnabledOXDNAFlag::OXRNA2;
}
}    // namespace LAMMPS_NS

#endif
#endif

