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

FixStyle(SHEAR_HISTORY,FixShearHistory)

#else

#ifndef LMP_FIX_SHEAR_HISTORY_H
#define LMP_FIX_SHEAR_HISTORY_H

#include "fix.h"
#include "my_page.h"

namespace LAMMPS_NS {

class FixShearHistory : public Fix {
  friend class Neighbor;
  friend class PairGranHookeHistory;

 public:
  FixShearHistory(class LAMMPS *, int, char **);
  ~FixShearHistory();
  int setmask();
  void init();
  void setup_pre_exchange();
  virtual void pre_exchange();
  void min_setup_pre_exchange();
  void min_pre_exchange();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 protected:
  int *npartner;                // # of touching partners of each atom
  tagint **partner;             // global atom IDs for the partners
  double (**shearpartner)[3];   // 3 shear values with the partner
  int maxtouch;                 // max # of touching partners for my atoms

  class Pair *pair;
  int *computeflag;             // computeflag in PairGranHookeHistory

  int pgsize,oneatom;           // copy of settings in Neighbor
  MyPage<tagint> *ipage;        // pages of partner atom IDs
  MyPage<double[3]> *dpage;     // pages of shear history with partners

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

UNDOCUMENTED

U: Too many touching neighbors - boost MAXTOUCH

A granular simulation has too many neighbors touching one atom.  The
MAXTOUCH parameter in fix_shear_history.cpp must be set larger and
LAMMPS must be re-built.

*/
