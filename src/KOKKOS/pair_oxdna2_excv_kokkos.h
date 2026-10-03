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
PairStyle(oxdna2/excv/kk,PairOxdna2ExcvKokkos<LMPDeviceType>);
PairStyle(oxdna2/excv/kk/device,PairOxdna2ExcvKokkos<LMPDeviceType>);
PairStyle(oxdna2/excv/kk/host,PairOxdna2ExcvKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA2_EXCV_KOKKOS_H
#define LMP_PAIR_OXDNA2_EXCV_KOKKOS_H

#include "pair_oxdna_excv_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxdna2ExcvKokkos : public PairOxdnaExcvKokkos<DeviceType> {
 public:
  PairOxdna2ExcvKokkos(class LAMMPS *);
  ~PairOxdna2ExcvKokkos() {}
};

template<class DeviceType>
PairOxdna2ExcvKokkos<DeviceType>::PairOxdna2ExcvKokkos(LAMMPS *lmp) : PairOxdnaExcvKokkos<DeviceType>(lmp)
{
   this->oxdnaflag = PairOxdnaExcvKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA2;
}
}    // namespace LAMMPS_NS

#endif
#endif
