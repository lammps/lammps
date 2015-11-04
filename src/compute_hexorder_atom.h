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

ComputeStyle(hexorder/atom,ComputeHexOrderAtom)

#else

#ifndef LMP_COMPUTE_HEXORDER_ATOM_H
#define LMP_COMPUTE_HEXORDER_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHexOrderAtom : public Compute {
 public:
  ComputeHexOrderAtom(class LAMMPS *, int, char **);
  ~ComputeHexOrderAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax,ncol;
  double cutsq;
  class NeighList *list;

  int *typelo,*typehi;
  double **q6array;

  void calc_q6(double, double, double&, double&);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute hexorder/atom requires a pair style be defined

Self-explantory.

E: Compute hexorder/atom cutoff is longer than pairwise cutoff

Cannot compute order parameter at distances longer than the pair cutoff,
since those atoms are not in the neighbor list.

W: More than one compute hexorder/atom

It is not efficient to use compute hexorder/atom more than once.

*/
