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

#ifndef COMPUTE_TEMP_RAMP_H
#define COMPUTE_TEMP_RAMP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRamp : public Compute {
 public:
  ComputeTempRamp(class LAMMPS *, int, char **);
  ~ComputeTempRamp();
  void init();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(double *);
  void restore_bias_all();
  double memory_usage();

 private:
  int coord_dim;
  double coord_lo,coord_hi;
  int v_dim;
  double v_lo,v_hi;
  int scaleflag,fix_dof;
  double tfactor,xscale,yscale,zscale;

  double vbias[3];    // stored velocity bias for one atom
  double **vbiasall;  // stored velocity bias for all atoms
  int maxbias;        // size of vbiasall array
  Compute *tbias;     // ptr to additional bias compute

  void recount();
};

}

#endif
