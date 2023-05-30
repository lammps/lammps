/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  ~FixNeighHistory() override;
  int setmask() override;
  void init() override;
  void setup_post_neighbor() override;
  void pre_exchange() override;
  void min_pre_exchange() override;
  void post_neighbor() override;
  void min_post_neighbor() override;
  void post_run() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;

  int pack_reverse_comm_size(int, int) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void write_restart(FILE *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

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

  // Shift by HISTBITS and check the first bit
  inline int histmask(int j) const { return j >> HISTBITS & 1; }

  enum { DEFAULT, NPARTNER, PERPARTNER };
};

}    // namespace LAMMPS_NS

#endif
#endif
