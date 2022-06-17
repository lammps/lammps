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

ComputeStyle(sna/grid/local,ComputeSNAGridLocal)

#else

#ifndef LMP_COMPUTE_SNA_GRID_LOCAL_H
#define LMP_COMPUTE_SNA_GRID_LOCAL_H

#include "compute_grid_local.h"

namespace LAMMPS_NS {

class ComputeSNAGridLocal : public ComputeGridLocal {
 public:
  ComputeSNAGridLocal(class LAMMPS *, int, char **);
  ~ComputeSNAGridLocal();
  void init();
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();

 private:
  int ncoeff;
  double **cutsq;
  class NeighList *list;
  double rcutfac;
  double *radelem;
  double *wjelem;
  int *map;    // map types to [0,nelements)
  int nelements, chemflag;
  int switchinnerflag;
  double *sinnerelem;
  double *dinnerelem;
  class SNA *snaptr;
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

W: More than one compute sna/grid/local

Self-explanatory.

*/
