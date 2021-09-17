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
ComputeStyle(snad/atom,ComputeSNADAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SNAD_ATOM_H
#define LMP_COMPUTE_SNAD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSNADAtom : public Compute {
 public:
  ComputeSNADAtom(class LAMMPS *, int, char **);
  ~ComputeSNADAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax;
  int ncoeff, nperdim, yoffset, zoffset;
  double **cutsq;
  class NeighList *list;
  double **snad;
  double rcutfac;
  double *radelem;
  double *wjelem;
  int *map;    // map types to [0,nelements)
  int nelements, chemflag;
  class SNA *snaptr;
  double cutmax;
  int quadraticflag;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute snad/atom requires a pair style be defined

Self-explanatory.

E: Compute sna/atom cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute snad/atom

Self-explanatory.

U: Compute snad/atom cutoff is longer than pairwise cutoff

Self-explanatory.

*/
