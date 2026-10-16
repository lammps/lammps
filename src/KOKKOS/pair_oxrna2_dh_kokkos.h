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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxrna2/dh/kk,PairOxrna2DhKokkos<LMPDeviceType>);
PairStyle(oxrna2/dh/kk/device,PairOxrna2DhKokkos<LMPDeviceType>);
PairStyle(oxrna2/dh/kk/host,PairOxrna2DhKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_DH_KOKKOS_H
#define LMP_PAIR_OXRNA2_DH_KOKKOS_H

#include "pair_oxdna2_dh_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxrna2DhKokkos : public PairOxdna2DhKokkos<DeviceType> {
 public:
  PairOxrna2DhKokkos(class LAMMPS *);
  ~PairOxrna2DhKokkos() {}
};

template<class DeviceType>
PairOxrna2DhKokkos<DeviceType>::PairOxrna2DhKokkos(LAMMPS *lmp) : PairOxdna2DhKokkos<DeviceType>(lmp)
{
   this->oxdnaflag = PairOxdna2DhKokkos<DeviceType>::EnabledOXDNAFlag::OXRNA2;
}
}    // namespace LAMMPS_NS

#endif
#endif

