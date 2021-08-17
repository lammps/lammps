/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(cnp/atom,ComputeCNPAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_CNP_ATOM_H
#define LMP_COMPUTE_CNP_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCNPAtom : public Compute {
 public:
  ComputeCNPAtom(class LAMMPS *, int, char **);
  ~ComputeCNPAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  //revise
  int nmax;
  double cutsq;
  class NeighList *list;
  int **nearest;
  int *nnearest;
  double *cnpv;
  //  int nmax;
  //  double cutsq;
  //  class NeighList *list;
  //  int **nearest;
  // int *nnearest;
  //  double *pattern;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute cnp/atom requires a pair style be defined

Self-explanatory.

E: Compute cnp/atom cutoff is longer than pairwise cutoff

Self-explanatory.

W: Compute cnp/atom cutoff may be too large to find ghost atom neighbors

The neighbor cutoff used may not encompass enough ghost atoms
to perform this operation correctly.

W: More than one compute cnp/atom defined

It is not efficient to use compute cnp/atom more than once.

W: Too many neighbors in CNP for %d atoms

More than the maximum # of neighbors was found multiple times. This
was unexpected.

W: Too many common neighbors in CNP %d times

More than the maximum # of neighbors was found multiple times. This
was unexpected.

*/
