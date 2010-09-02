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

FixStyle(nve/eff,FixNVEEff)

#else

#ifndef LMP_FIX_NVE_EFF_H
#define LMP_FIX_NVE_EFF_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEEff : public FixNVE {
 public:
  FixNVEEff(class LAMMPS *, int, char **);
  void initial_integrate(int);
  void final_integrate();
};

}

#endif
#endif
