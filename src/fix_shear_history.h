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

#ifdef FIX_CLASS

FixStyle(SHEAR_HISTORY,FixShearHistory)

#else

#ifndef LMP_FIX_SHEAR_HISTORY_H
#define LMP_FIX_SHEAR_HISTORY_H

#include "fix.h"
#include "my_page.h"

namespace LAMMPS_NS {

class FixShearHistory : public Fix {
  //friend class Neighbor;
  //friend class PairGranHookeHistory;
  friend class PairLineGranHookeHistory;
  friend class PairTriGranHookeHistory;

 public:
  int nlocal_neigh;             // nlocal at last time neigh list was built
  int nall_neigh;               // ditto for nlocal+nghost
  int *npartner;                // # of touching partners of each atom
  tagint **partner;             // global atom IDs for the partners
  double **shearpartner;        // shear values with the partner
  class Pair *pair;             // ptr to pair style that uses shear history

  FixShearHistory(class LAMMPS *, int, char **);
  ~FixShearHistory();
  int setmask();
  void init();
  virtual void pre_exchange();
  void min_pre_exchange();
  void post_run();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);

  int pack_reverse_comm_size(int, int);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 protected:
  int newton_pair;              // same as force setting
  int dnum,dnumbytes;           // dnum = # of shear history values
  int onesided;                 // 1 for line/tri history, else 0

  int maxtouch;                 // max # of touching partners for my atoms
  int commflag;                 // mode of reverse comm to get ghost info

  int pgsize,oneatom;           // copy of settings in Neighbor
  MyPage<tagint> *ipage;        // pages of partner atom IDs
  MyPage<double> *dpage;        // pages of shear history with partners

  void pre_exchange_onesided();
  void pre_exchange_newton();
  void pre_exchange_no_newton();
  void allocate_pages();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair style granular with history requires atoms have IDs

Atoms in the simulation do not have IDs, so history effects
cannot be tracked by the granular pair potential.

E: Shear history overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

*/
