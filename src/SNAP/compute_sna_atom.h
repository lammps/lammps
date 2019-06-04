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

ComputeStyle(sna/atom,ComputeSNAAtom)

#else

#ifndef LMP_COMPUTE_SNA_ATOM_H
#define LMP_COMPUTE_SNA_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSNAAtom : public Compute {
 public:
  ComputeSNAAtom(class LAMMPS *, int, char **);
  ~ComputeSNAAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax, njmax, diagonalstyle;
  int ncoeff;
  double **cutsq;
  class NeighList *list;
  double **sna;
  double rcutfac;
  double *radelem;
  double *wjelem;
  class SNA* snaptr;
  double cutmax;
  int quadraticflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute sna/atom requires a pair style be defined

Self-explanatory.

E: Compute sna/atom cutoff is longer than pairwise cutoff

Self-explanatory.

W: More than one compute sna/atom

Self-explanatory.

*/
