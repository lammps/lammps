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
PairStyle(oxdna3/hbond/kk,PairOxdna3HbondKokkos<LMPDeviceType>);
PairStyle(oxdna3/hbond/kk/device,PairOxdna3HbondKokkos<LMPDeviceType>);
PairStyle(oxdna3/hbond/kk/host,PairOxdna3HbondKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_HBOND_KOKKOS_H
#define LMP_PAIR_OXDNA3_HBOND_KOKKOS_H

#include "pair_oxdna_hbond_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxdna3HbondKokkos : public PairOxdnaHbondKokkos<DeviceType> {
 public:
  PairOxdna3HbondKokkos(class LAMMPS *);
  ~PairOxdna3HbondKokkos() {}
  void coeff(int, char **) override;
};

/* ----------------------------------------------------------------------
   IMPORTANT NOTE ! We entirely code duplicate the sequence-specific alpha_hb
   setup between PairOxdna3HbondKokkos and PairOxdna3Hbond. So any edits made
   in one need to manually be made to the other !
   The vanilla version is in: src/CG-DNA/pair_oxdna3_hbond.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna3HbondKokkos<DeviceType>::PairOxdna3HbondKokkos(LAMMPS *lmp) : PairOxdnaHbondKokkos<DeviceType>(lmp)
{
  this->oxdnaflag = PairOxdnaHbondKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA3;

  // sequence-specific base-pairing strength
  // A:0 C:1 G:2 T:3, 5'- [i][j] -3'

  this->alpha_hb[0][0] = 1.00000;
  this->alpha_hb[0][1] = 1.00000;
  this->alpha_hb[0][2] = 1.00000;
  this->alpha_hb[0][3] = 0.6493620379646540;

  this->alpha_hb[1][0] = 1.00000;
  this->alpha_hb[1][1] = 1.00000;
  this->alpha_hb[1][2] = 1.1999420813642658;
  this->alpha_hb[1][3] = 1.00000;

  this->alpha_hb[2][0] = 1.00000;
  this->alpha_hb[2][1] = 1.1999420813642658;
  this->alpha_hb[2][2] = 1.00000;
  this->alpha_hb[2][3] = 1.00000;

  this->alpha_hb[3][0] = 0.6493620379646540;
  this->alpha_hb[3][1] = 1.00000;
  this->alpha_hb[3][2] = 1.00000;
  this->alpha_hb[3][3] = 1.00000;
}

}    // namespace LAMMPS_NS

#endif
#endif
