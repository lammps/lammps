/* -*- c++ -*- ----------------------------------------------------------
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

ComputeStyle(cac/quad/count,ComputeCACQuadCount)

#else

#ifndef LMP_COMPUTE_CAC_QUAD_COUNT_H
#define LMP_COMPUTE_CAC_QUAD_COUNT_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCACQuadCount : public Compute {
 public:
	 ComputeCACQuadCount(class LAMMPS *, int, char **);
  ~ComputeCACQuadCount();
  void init();
  void compute_peratom();
  double memory_usage();
  void compute_surface_depths(double &x, double &y, double &z,
	  int &xb, int &yb, int &zb, int flag);
 private:
  int nmax;
  double *quad_count;
  double ***current_nodal_positions;
  int current_element_type, current_poly_count;
  int *current_element_scale;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: weight counts don't match nlocal

Internal debug check in case assumptions break down for a given simulation input. Contact Author for help.

E: compute cac/quad/count requires a CAC atom style

Self-explanatory

W: More than one compute cac/quad/count

It is redundant to define more than one of these computes given their use in tandem with fix balance to load balance models

*/
