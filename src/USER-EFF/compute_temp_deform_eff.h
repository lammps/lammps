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

ComputeStyle(temp/deform/eff,ComputeTempDeformEff)

#else

#ifndef LMP_COMPUTE_TEMP_DEFORM_EFF_H
#define LMP_COMPUTE_TEMP_DEFORM_EFF_H

#include "compute_temp_deform.h"

namespace LAMMPS_NS {

class ComputeTempDeformEff : public ComputeTempDeform {
 public:
  ComputeTempDeformEff(class LAMMPS *, int, char **);
  ~ComputeTempDeformEff() {}
  double compute_scalar();
  void compute_vector();

 private:
  void dof_compute();
};

}

#endif
#endif
