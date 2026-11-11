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
PairStyle(vashishta/table/kk,PairVashishtaTableKokkos<LMPDeviceType>);
PairStyle(vashishta/table/kk/device,PairVashishtaTableKokkos<LMPDeviceType>);
PairStyle(vashishta/table/kk/host,PairVashishtaTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_VASHISHTA_TABLE_KOKKOS_H
#define LMP_PAIR_VASHISHTA_TABLE_KOKKOS_H

#include "pair_vashishta_kokkos.h"

namespace LAMMPS_NS {

// same kernels as vashishta/kk; the two-body term switches to linear
// interpolation in rsq beyond the inner cutoff through "tabflag"

template<class DeviceType>
class PairVashishtaTableKokkos : public PairVashishtaKokkos<DeviceType> {
 public:
  PairVashishtaTableKokkos(class LAMMPS *);
  void settings(int, char **) override;
  double memory_usage() override;

 protected:
  int ntable;
  DAT::tdual_kkfloat_3d k_forceTable, k_potentialTable;

  void setup_params() override;
};

}

#endif
#endif
