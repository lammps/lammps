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
PairStyle(oxdna3/coaxstk/kk,PairOxdna3CoaxstkKokkos<LMPDeviceType>);
PairStyle(oxdna3/coaxstk/kk/device,PairOxdna3CoaxstkKokkos<LMPDeviceType>);
PairStyle(oxdna3/coaxstk/kk/host,PairOxdna3CoaxstkKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_COAXSTK_KOKKOS_H
#define LMP_PAIR_OXDNA3_COAXSTK_KOKKOS_H

#include "pair_oxdna2_coaxstk_kokkos.h"
#include "pair_oxdna3_coaxstk.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairOxdna3CoaxstkKokkos : public PairOxdna2CoaxstkKokkos<DeviceType> {
 public:
  PairOxdna3CoaxstkKokkos(class LAMMPS *);
  ~PairOxdna3CoaxstkKokkos() {}
};

/* ----------------------------------------------------------------------
   IMPORTANT NOTE ! We entirely code duplicate the sequence-specific eta_cxst
   setup between PairOxdna3CoaxstkKokkos and PairOxdna3Coaxstk. So any edits
   made in one need to manually be made to the other !
   The vanilla version is in: src/CG-DNA/pair_oxdna3_coaxstk.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna3CoaxstkKokkos<DeviceType>::PairOxdna3CoaxstkKokkos(LAMMPS *lmp) : PairOxdna2CoaxstkKokkos<DeviceType>(lmp)
{
  this->oxdnaflag = PairOxdna2CoaxstkKokkos<DeviceType>::EnabledOXDNAFlag::OXDNA3;

  // sequence-specific coaxial stacking strength
  // A:0 C:1 G:2 T:3, 3'- [i] X [j] -5'

  this->eta_cxst[0][0] = 1.1217958408368172;
  this->eta_cxst[1][0] = 1.0712851690057155;
  this->eta_cxst[2][0] = 1.1161603311902566;
  this->eta_cxst[3][0] = 1.0052361315065244;

  this->eta_cxst[0][1] = 1.1217958408368172;
  this->eta_cxst[1][1] = 0.7892685731520542;
  this->eta_cxst[2][1] = 1.1022201982984874;
  this->eta_cxst[3][1] = 0.8658975520778347;

  this->eta_cxst[0][2] = 1.1217958408368172;
  this->eta_cxst[1][2] = 0.9896542231533637;
  this->eta_cxst[2][2] = 1.1088392608169480;
  this->eta_cxst[3][2] = 1.1217958408368172;

  this->eta_cxst[0][3] = 0.9300223683636719;
  this->eta_cxst[1][3] = 0.7694592613578328;
  this->eta_cxst[2][3] = 1.0007533199170144;
  this->eta_cxst[3][3] = 0.8593983791552220;
}

}    // namespace LAMMPS_NS

#endif
#endif
