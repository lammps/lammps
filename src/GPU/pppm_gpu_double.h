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

KSpaceStyle(pppm/gpu/double,PPPMGPUDouble)

#else

#ifndef LMP_PPPM_GPU_DOUBLE_H
#define LMP_PPPM_GPU_DOUBLE_H

#include "pppm_gpu.h"
#include "lmptype.h"

namespace LAMMPS_NS {

class PPPMGPUDouble : public PPPMGPU<double> {
 public:
  PPPMGPUDouble(class LAMMPS *, int, char **);
  ~PPPMGPUDouble();
  void init();
  void compute(int, int);
  double memory_usage();

 protected:
};

}

#endif
#endif
