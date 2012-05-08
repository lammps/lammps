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

#ifdef FIX_CLASS

FixStyle(reax/c/bonds,FixReaxCBonds)

#else

#ifndef LMP_FIX_REAXC_BONDS_H
#define LMP_FIX_REAXC_BONDS_H

#include "stdio.h"
#include "fix.h"
#include "pair_reax_c.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"
#include "pointers.h"

#define MAXBOND 24

namespace LAMMPS_NS {

class FixReaxCBonds : public Fix {
 public:
  FixReaxCBonds(class LAMMPS *, int, char **);
  ~FixReaxCBonds();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
  int me, nprocs, nmax, ntypes, maxsize;
  int nrepeat, irepeat, repeat, nfreq;
  int *numneigh, **neighid, **tmpid;
  double *sbo, *nlp, *avq, **abo, **tmpabo;
  FILE *fp;

  void allocate();
  void Output_ReaxC_Bonds(bigint, FILE *);
  void GatherBond(reax_system*, reax_list*);
  void FindBond(reax_system*, reax_list*, int &);
  void PassBuffer(reax_system*, double *, int &);
  void RecvBuffer(reax_system*, double *, int, int, int, int);
  int nint(const double &);
  double memory_usage();

  bigint nvalid, nextvalid();
  reax_system *system;
  reax_list *lists;
  class PairReaxC *reaxc;
  class NeighList *list;

};
}

#endif
#endif
