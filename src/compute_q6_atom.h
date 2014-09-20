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

ComputeStyle(q6/atom,ComputeQ6Atom)

#else

#ifndef LMP_COMPUTE_Q6_ATOM_H
#define LMP_COMPUTE_Q6_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeQ6Atom : public Compute {
 public:
  ComputeQ6Atom(class LAMMPS *, int, char **);
  ~ComputeQ6Atom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double q6cutoff;
  class NeighList *list;
  double *q6;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: More than one compute q6/atom

It is not efficient to use compute q6/atom more than once.

*/
