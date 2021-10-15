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

ComputeStyle(phase/atom,ComputePhaseAtom)

#else

#ifndef LMP_COMPUTE_PHASE_ATOM_H
#define LMP_COMPUTE_PHASE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePhaseAtom : public Compute {
 public:
  ComputePhaseAtom(class LAMMPS *, int, char **);
  virtual ~ComputePhaseAtom();
  virtual void init();
  void init_list(int, class NeighList *);
  virtual void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 protected:
  int nmax;
  double cutoff,cutsq,sphere_vol;
  class NeighList *list;

  double **phase;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute phase/atom compute ID

UNDOCUMENTED

E: Compute phase/atom compute ID is not orientorder/atom

UNDOCUMENTED

E: Compute phase/atom threshold not between -1 and 1

UNDOCUMENTED

E: Invalid cstyle in compute phase/atom

UNDOCUMENTED

E: Compute phase/atom requires components option in compute orientorder/atom

UNDOCUMENTED

E: Compute phase/atom requires a pair style be defined

Self-explanatory.

E: Compute phase/atom cutoff is longer than pairwise cutoff

Cannot compute phase at distances longer than the pair cutoff,
since those atoms are not in the neighbor list.

*/
