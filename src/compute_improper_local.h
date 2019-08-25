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

ComputeStyle(improper/local,ComputeImproperLocal)

#else

#ifndef LMP_COMPUTE_IMPROPER_LOCAL_H
#define LMP_COMPUTE_IMPROPER_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeImproperLocal : public Compute {
 public:
  ComputeImproperLocal(class LAMMPS *, int, char **);
  ~ComputeImproperLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues,cflag;
  int ncount;

  int nmax;
  double *vlocal;
  double **alocal;

  int compute_impropers(int);
  void reallocate(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute improper/local used when impropers are not allowed

The atom style does not support impropers.

E: Invalid keyword in compute improper/local command

Self-explanatory.

E: No improper style is defined for compute improper/local

Self-explanatory.

*/
