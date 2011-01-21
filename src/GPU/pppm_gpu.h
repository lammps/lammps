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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/gpu,PPPM_GPU)

#else

#ifndef LMP_PPPM_GPU_H
#define LMP_PPPM_GPU_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPM_GPU : public PPPM {
 public:
  PPPM_GPU(class LAMMPS *, int, char **);
  ~PPPM_GPU();
  void init();
  void setup();
  void compute(int, int);

 protected:
};

}

#endif
#endif
