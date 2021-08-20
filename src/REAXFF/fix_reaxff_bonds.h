// clang-format off
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

#ifdef FIX_CLASS
// clang-format off
FixStyle(reaxff/bonds,FixReaxFFBonds);
FixStyle(reax/c/bonds,FixReaxFFBonds);
// clang-format on
#else

#ifndef LMP_FIX_REAXC_BONDS_H
#define LMP_FIX_REAXC_BONDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReaxFFBonds : public Fix {
 public:
  FixReaxFFBonds(class LAMMPS *, int, char **);
  virtual ~FixReaxFFBonds();
  int setmask();
  virtual void init();
  void setup(int);
  void end_of_step();

 protected:
  int me, nprocs, nmax, ntypes, maxsize;
  int *numneigh;
  tagint **neighid;
  double **abo;
  FILE *fp;

  void allocate();
  void destroy();
  virtual void Output_ReaxFF_Bonds();
  int FindBond();
  void PassBuffer(double *, int &);
  void RecvBuffer(double *, int, int, int, int);
  int nint(const double &);
  virtual double memory_usage();

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxFF *reaxff;
  class NeighList *list;
};
}

#endif
#endif
