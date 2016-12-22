/* ----------------------------------------------------------------------
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

FixStyle(dpd/energy/kk,FixDPDenergyKokkos<LMPDeviceType>)
FixStyle(dpd/energy/kk/device,FixDPDenergyKokkos<LMPDeviceType>)
FixStyle(dpd/energy/kk/host,FixDPDenergyKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_DPDE_H
#define LMP_FIX_DPDE_H

#include "fix_dpd_energy.h"

namespace LAMMPS_NS {

class FixDPDenergyKokkos : public FixDPDEnergy {
 public:
  FixDPDenergyKokkos(class LAMMPS *, int, char **);
  virtual ~FixDPDenergyKokkos() {}
  virtual void initial_integrate(int);
  virtual void final_integrate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
