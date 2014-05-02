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

ComputeStyle(temp/region,ComputeTempRegion)

#else

#ifndef LMP_COMPUTE_TEMP_REGION_H
#define LMP_COMPUTE_TEMP_REGION_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRegion : public Compute {
 public:
  ComputeTempRegion(class LAMMPS *, int, char **);
  virtual ~ComputeTempRegion();
  void init();
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();

  void dof_remove_pre();
  int dof_remove(int);
  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 protected:
  int iregion;
  char *idregion;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for compute temp/region does not exist

Self-explanatory.

*/
