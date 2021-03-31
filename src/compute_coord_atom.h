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

ComputeStyle(coord/atom,ComputeCoordAtom)

#else

#ifndef LMP_COMPUTE_COORD_ATOM_H
#define LMP_COMPUTE_COORD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCoordAtom : public Compute {
 public:
  ComputeCoordAtom(class LAMMPS *, int, char **);
  virtual ~ComputeCoordAtom();
  virtual void init();
  void init_list(int, class NeighList *);
  virtual void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();
  enum {NONE,CUTOFF,ORIENT};

 protected:
  int nmax,ncol;
  double cutsq;
  class NeighList *list;

  int *typelo,*typehi;
  double *cvec;
  double **carray;

  char *group2;
  int jgroup,jgroupbit;

  class ComputeOrientOrderAtom *c_orientorder;
  char *id_orientorder;
  double threshold;
  double **normv;
  int cstyle,nqlist,l;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute coord/atom compute ID

UNDOCUMENTED

E: Compute coord/atom compute ID is not orientorder/atom

UNDOCUMENTED

E: Compute coord/atom threshold not between -1 and 1

UNDOCUMENTED

E: Invalid cstyle in compute coord/atom

UNDOCUMENTED

E: Compute coord/atom requires components option in compute orientorder/atom

UNDOCUMENTED

E: Compute coord/atom requires a pair style be defined

Self-explanatory.

E: Compute coord/atom cutoff is longer than pairwise cutoff

Cannot compute coordination at distances longer than the pair cutoff,
since those atoms are not in the neighbor list.

*/
