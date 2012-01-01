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

ComputeStyle(slice,ComputeSlice)

#else

#ifndef LMP_COMPUTE_SLICE_H
#define LMP_COMPUTE_SLICE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSlice : public Compute {
 public:
  ComputeSlice(class LAMMPS *, int, char **);
  virtual ~ComputeSlice();
  void init();
  void compute_vector();
  void compute_array();

 private:
  int me;
  int nstart,nstop,nskip,nvalues;
  int *which,*argindex,*value2index;
  char **ids;

  void extract_one(int, double *, int);
};

}

#endif
#endif
