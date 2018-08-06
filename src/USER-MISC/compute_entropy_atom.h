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

ComputeStyle(entropy/atom,ComputeEntropyAtom)

#else

#ifndef COMPUTE_ENTROPY_ATOM_H
#define COMPUTE_ENTROPY_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEntropyAtom : public Compute {
 public:
  ComputeEntropyAtom(class LAMMPS *, int, char **);
  ~ComputeEntropyAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax,maxneigh, nbin;
  class NeighList *list;
  double *pair_entropy, *pair_entropy_avg;
  double sigma, cutoff, cutoff2;
  double cutsq, cutsq2;
  double deltar;
  int deltabin;
  int avg_flag;
  int local_flag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute entropy/atom requires a pair style be defined

This is because the computation of the pair entropy values
uses a pairwise neighbor list.

W: More than one compute entropy/atom

It is not efficient to use compute entropy/atom more than once.

*/
