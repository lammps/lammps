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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/eff,ComputeTempEff)

#else

#ifndef LMP_COMPUTE_TEMP_EFF_H
#define LMP_COMPUTE_TEMP_EFF_H

#include "compute.h"

namespace LAMMPS_NS {
	
class ComputeTempEff : public Compute {
 public:
  ComputeTempEff(class LAMMPS *, int, char **);
  virtual ~ComputeTempEff();
  void init();
  double compute_scalar();
  void compute_vector();
  
 private:
  int fix_dof;
  double tfactor;
  double *inertia;

  void dof_compute();
};
	
}

#endif
#endif
