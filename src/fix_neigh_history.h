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
FixStyle(NEIGH_HISTORY,FixNeighHistory);
// clang-format on
#else

#ifndef LMP_FIX_NEIGH_HISTORY_H
#define LMP_FIX_NEIGH_HISTORY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNeighHistory : public Fix {
 public:
  int nlocal_neigh;       // nlocal at last time neigh list was built
  int nall_neigh;         // ditto for nlocal+nghost
  int use_bit_flag;       // flag whether this fix uses the extra bit in the nlist
  int **firstflag;        // ptr to each atom's neighbor flsg
  double **firstvalue;    // ptr to each atom's values
  class Pair *pair;       // ptr to pair style that uses neighbor history

  FixNeighHistory(class LAMMPS *, int, char **);
  ~FixNeighHistory();
  int setmask();
  void init();
  void setup_post_neighbor();
  void pre_exchange();
  void min_pre_exchange();
  virtual void post_neighbor();
  void min_post_neighbor();
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
  void write_restart(FILE *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 protected:
  int newton_pair;        // same as force setting
  int dnum, dnumbytes;    // dnum = # of values per neighbor
  int onesided;           // 1 for line/tri history, else 0

  int maxatom;     // max size of firstflag and firstvalue
  int commflag;    // mode of reverse comm to get ghost info
  double *zeroes;

  // per-atom data structures
  // partners = flagged neighbors of an atom

  int *npartner;            // # of partners of each atom
  tagint **partner;         // global atom IDs for the partners
  double **valuepartner;    // values for the partners
  int maxpartner;           // max # of partners for any of my atoms

  // per-atom data structs pointed to by partner & valuepartner

  int pgsize, oneatom;           // copy of settings in Neighbor
  MyPage<tagint> *ipage_atom;    // pages of partner atom IDs
  MyPage<double> *dpage_atom;    // pages of partner values

  // per-neighbor data structs pointed to by firstflag & firstvalue

  MyPage<int> *ipage_neigh;       // pages of local atom indices
  MyPage<double> *dpage_neigh;    // pages of partner values

  virtual void pre_exchange_onesided();
  virtual void pre_exchange_newton();
  virtual void pre_exchange_no_newton();
  void allocate_pages();

  inline int sbmask(int j) const { return j >> SBBITS & 3; }
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Neighbor history requires atoms have IDs

UNDOCUMENTED

E: Neighbor history overflow, boost neigh_modify one

UNDOCUMENTED

E: Unsupported comm mode in neighbor history

UNDOCUMENTED

U: Pair style granular with history requires atoms have IDs

Atoms in the simulation do not have IDs, so history effects
cannot be tracked by the granular pair potential.

U: Shear history overflow, boost neigh_modify one

There are too many neighbors of a single atom.  Use the neigh_modify
command to increase the max number of neighbors allowed for one atom.
You may also want to boost the page size.

*/
