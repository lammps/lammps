/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(dpd/energy/kk,FixDPDenergyKokkos<LMPDeviceType>);
FixStyle(dpd/energy/kk/device,FixDPDenergyKokkos<LMPDeviceType>);
FixStyle(dpd/energy/kk/host,FixDPDenergyKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DPDE_KOKKOS_H
#define LMP_FIX_DPDE_KOKKOS_H

#include "fix_dpd_energy.h"
#include "pair_dpd_fdt_energy_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixDPDenergyKokkos : public FixDPDenergy {
 public:
  FixDPDenergyKokkos(class LAMMPS *, int, char **);
  virtual ~FixDPDenergyKokkos() {}
  virtual void initial_integrate(int);
  virtual void final_integrate();

  void take_half_step();
 protected:
  PairDPDfdtEnergyKokkos<DeviceType>* pairDPDEKK;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must use pair_style dpd/fdt/energy/kk with fix dpd/energy/kk

Self-explanatory.

*/
