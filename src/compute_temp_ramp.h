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

ComputeStyle(temp/ramp,ComputeTempRamp)

#else

#ifndef LMP_COMPUTE_TEMP_RAMP_H
#define LMP_COMPUTE_TEMP_RAMP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRamp : public Compute {
 public:
  ComputeTempRamp(class LAMMPS *, int, char **);
  ~ComputeTempRamp();
  void init() {}
  void setup();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 private:
  int coord_dim;
  double coord_lo,coord_hi;
  int v_dim;
  double v_lo,v_hi;
  int scaleflag,fix_dof;
  double tfactor,xscale,yscale,zscale;

  void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Use of compute temp/ramp with undefined lattice

Must use lattice command with compute temp/ramp command if units
option is set to lattice.

*/
