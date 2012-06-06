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

ComputeStyle(cluster/atom,ComputeClusterAtom)

#else

#ifndef LMP_COMPUTE_CLUSTER_ATOM_H
#define LMP_COMPUTE_CLUSTER_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeClusterAtom : public Compute {
 public:
  ComputeClusterAtom(class LAMMPS *, int, char **);
  ~ComputeClusterAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  double memory_usage();

 private:
  int nmax;
  double cutsq;
  class NeighList *list;
  double *clusterID;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use compute cluster/atom unless atoms have IDs

Atom IDs are used to identify clusters.

E: Compute cluster/atom requires a pair style be defined

This is so that the pair style defines a cutoff distance which
is used to find clusters.

E: Compute cluster/atom cutoff is longer than pairwise cutoff

Cannot identify clusters beyond cutoff.

W: More than one compute cluster/atom

It is not efficient to use compute cluster/atom  more than once.

*/
