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
PairStyle(oxdna3/excv/kk,PairOxdna3ExcvKokkos<LMPDeviceType>);
PairStyle(oxdna3/excv/kk/device,PairOxdna3ExcvKokkos<LMPDeviceType>);
PairStyle(oxdna3/excv/kk/host,PairOxdna3ExcvKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_EXCV_KOKKOS_H
#define LMP_PAIR_OXDNA3_EXCV_KOKKOS_H

#include "pair_oxdna_excv_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxdna3ExcvKokkos : public PairOxdnaExcvKokkos<DeviceType> {
 public:
  PairOxdna3ExcvKokkos(class LAMMPS *);
  ~PairOxdna3ExcvKokkos() {}
  void coeff(int, char **) override;
};

template<class DeviceType>
PairOxdna3ExcvKokkos<DeviceType>::PairOxdna3ExcvKokkos(LAMMPS *lmp) : PairOxdnaExcvKokkos<DeviceType>(lmp)
{
    this->oxdnaflag = PairOxdnaExcvKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA3;
}

}    // namespace LAMMPS_NS

#endif
#endif
