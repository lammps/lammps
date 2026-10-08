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
PairStyle(oxrna2/hbond/kk,PairOxrna2HbondKokkos<LMPDeviceType>);
PairStyle(oxrna2/hbond/kk/device,PairOxrna2HbondKokkos<LMPDeviceType>);
PairStyle(oxrna2/hbond/kk/host,PairOxrna2HbondKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_HBOND_KOKKOS_H
#define LMP_PAIR_OXRNA2_HBOND_KOKKOS_H

#include "pair_oxdna_hbond_kokkos.h"
#include "pair_oxrna2_hbond.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxrna2HbondKokkos : public PairOxdnaHbondKokkos<DeviceType> {
 public:
  PairOxrna2HbondKokkos(class LAMMPS *);
  ~PairOxrna2HbondKokkos() {}
};

/* ----------------------------------------------------------------------
   IMPORTANT NOTE ! We entirely code duplicate the sequence-specific alpha_hb
   setup between PairOxrna2HbondKokkos and PairOxrna2Hbond. So any edits made
   in one need to manually be made to the other !
   The vanilla version is in: src/CG-DNA/pair_oxrna2_hbond.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
PairOxrna2HbondKokkos<DeviceType>::PairOxrna2HbondKokkos(LAMMPS *lmp) :
  PairOxdnaHbondKokkos<DeviceType>(lmp)
{
  // sequence-specific base-pairing strength
  // A:0 C:1 G:2 U:3, 5'- [i][j] -3'

  this->alpha_hb[0][0] = 1.00000;
  this->alpha_hb[0][1] = 1.00000;
  this->alpha_hb[0][2] = 1.00000;
  this->alpha_hb[0][3] = 0.94253;

  this->alpha_hb[1][0] = 1.00000;
  this->alpha_hb[1][1] = 1.00000;
  this->alpha_hb[1][2] = 1.22288;
  this->alpha_hb[1][3] = 1.00000;

  this->alpha_hb[2][0] = 1.00000;
  this->alpha_hb[2][1] = 1.22288;
  this->alpha_hb[2][2] = 1.00000;
  this->alpha_hb[2][3] = 0.58655;

  this->alpha_hb[3][0] = 0.94253;
  this->alpha_hb[3][1] = 1.00000;
  this->alpha_hb[3][2] = 0.58655;
  this->alpha_hb[3][3] = 1.00000;
}

}    // namespace LAMMPS_NS

#endif
#endif
