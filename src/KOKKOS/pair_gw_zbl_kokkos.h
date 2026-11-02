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
PairStyle(gw/zbl/kk,PairGWZBLKokkos<LMPDeviceType>);
PairStyle(gw/zbl/kk/device,PairGWZBLKokkos<LMPDeviceType>);
PairStyle(gw/zbl/kk/host,PairGWZBLKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GW_ZBL_KOKKOS_H
#define LMP_PAIR_GW_ZBL_KOKKOS_H

#include "pair_gw_kokkos.h"
#include "pair_gw_zbl.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairGWZBLKokkos : public PairGWKokkosBase<DeviceType,PairGWZBL> {
 public:
  PairGWZBLKokkos(class LAMMPS *lmp) : PairGWKokkosBase<DeviceType,PairGWZBL>(lmp)
  {
    this->zblflag = 1;
    this->zbl_a0 = this->global_a_0;
    this->zbl_epsilon0 = this->global_epsilon_0;
    this->zbl_e = this->global_e;
  }
};

}

#endif
#endif
