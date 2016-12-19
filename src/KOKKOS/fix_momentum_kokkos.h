/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(momentum/kk,FixMomentumKokkos<LMPDeviceType>)
FixStyle(momentum/kk/device,FixMomentumKokkos<LMPDeviceType>)
FixStyle(momentum/kk/host,FixMomentumKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_MOMENTUM_KOKKOS_H
#define LMP_FIX_MOMENTUM_KOKKOS_H

#include "fix_momentum.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixMomentumKokkos : public FixMomentum {
 public:
  typedef ArrayTypes<DeviceType> AT;

  FixMomentum(class LAMMPS *, int, char **);
  void init();
  void end_of_step();

 private:
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
