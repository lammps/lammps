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

FixStyle(SHEAR_HISTORY/omp,FixShearHistoryOMP)

#else

#ifndef LMP_FIX_SHEAR_HISTORY_OMP_H
#define LMP_FIX_SHEAR_HISTORY_OMP_H

#include "fix_shear_history.h"

namespace LAMMPS_NS {

class FixShearHistoryOMP : public FixShearHistory {

 public:
  FixShearHistoryOMP(class LAMMPS *lmp, int narg, char **argv)
    : FixShearHistory(lmp,narg,argv) {};
  virtual void pre_exchange();
};

}

#endif
#endif
