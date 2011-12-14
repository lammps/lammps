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

FixStyle(PERI_NEIGH_OMP,FixPeriNeighOMP)

#else

#ifndef LMP_FIX_PERI_NEIGH_OMP_H
#define LMP_FIX_PERI_NEIGH_OMP_H

#include "fix_peri_neigh.h"

namespace LAMMPS_NS {

class FixPeriNeighOMP : public FixPeriNeigh {

 public:
  FixPeriNeighOMP(class LAMMPS *lmp, int narg, char **argv) : 
    FixPeriNeigh(lmp,narg,argv) {};
  virtual void init();
};

}

#endif
#endif
