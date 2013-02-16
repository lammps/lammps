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

ComputeStyle(voronoi/atom,ComputeVoronoi)

#else

#ifndef LMP_COMPUTE_VORONOI_H
#define LMP_COMPUTE_VORONOI_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVoronoi : public Compute {
 public:
  ComputeVoronoi(class LAMMPS *, int, char **);
  ~ComputeVoronoi();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double **voro;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute voronoi/atom not allowed for triclinic boxes

This is a current restriction of this command.

W: More than one compute voronoi/atom command

It is not efficient to use compute voronoi/atom more than once.

*/
