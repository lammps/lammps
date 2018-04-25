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

ComputeStyle(aggregate/atom,ComputeAggregateAtom)

#else

#ifndef LMP_COMPUTE_AGGREGATE_ATOM_H
#define LMP_COMPUTE_AGGREGATE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAggregateAtom : public Compute {
 public:
  ComputeAggregateAtom(class LAMMPS *, int, char **);
  ~ComputeAggregateAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax,commflag;
  double cutsq;
  class NeighList *list;
  double *aggregateID;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute aggregate/atom used when bonds are not allowed

UNDOCUMENTED

E: Cannot use compute aggregate/atom unless atoms have IDs

Atom IDs are used to identify aggregates.

E: Compute aggregate/atom requires a bond style to be defined

This is so that a bond list is generated which is used to find aggregates.

E: Compute cluster/atom requires a pair style to be defined

UNDOCUMENTED

E: Compute cluster/atom cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute aggregate/atom

It is not efficient to use compute aggregate/atom  more than once.

*/
