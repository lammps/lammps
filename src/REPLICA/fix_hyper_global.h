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
FixStyle(hyper/global,FixHyperGlobal);
// clang-format on
#else

#ifndef LMP_FIX_HYPER_GLOBAL_H
#define LMP_FIX_HYPER_GLOBAL_H

#include "fix_hyper.h"

namespace LAMMPS_NS {

class FixHyperGlobal : public FixHyper {
 public:
  FixHyperGlobal(class LAMMPS *, int, char **);
  ~FixHyperGlobal();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup_pre_neighbor();
  void setup_pre_reverse(int, int);
  void pre_neighbor();
  void pre_reverse(int, int);
  double compute_scalar();
  double compute_vector(int);
  double query(int);

  double memory_usage();

  // extra methods visible to callers

  void init_hyper();
  void build_bond_list(int);

 private:
  int me;
  double cutbond, qfactor, vmax, tequil;

  int firstflag, bcastflag, owner, nevent, nevent_atom;
  double cutbondsq, beta, dt, t_hyper, invqfactorsq;
  double outvec[5];     // same as VECLEN in *.cpp
  double maxbondlen;    // max length of any bond
  double maxdriftsq;    // max distance any atom drifts from original pos
  int nobias;           // # of steps when bias = 0, b/c bond too long
  int negstrain;        // # of steps when biased bond has negative strain
  bigint groupatoms;    // # of atoms in fix group

  class NeighList *list;

  // list of my owned bonds
  // persists on a proc from one event until the next

  int maxbond;    // allocated size of blist

  struct OneBond {     // single IJ bond, atom I is owner
    int i, j;          // current local indices of 2 bond atoms
    int iold, jold;    // local indices when bonds were formed
    double r0;         // relaxed bond length
  };

  OneBond *blist;    // list of owned bonds
  int nblocal;       // # of owned bonds

  // coords and IDs of owned+ghost atoms when bonds were formed
  // persists on a proc from one event until the next

  int nlocal_old;    // nlocal for old atoms
  int nall_old;      // nlocal+nghost for old atoms
  int maxold;        // allocated size of old atoms

  double **xold;     // coords of atoms when bonds were formed
  tagint *tagold;    // IDs of atoms when bonds were forme
  int *old2now;      // o2n[i] = current local index of old atom I

  // MPI data struct for finding bond with max strain via Allreduce

  struct Two {
    double value;
    int proc;
  };
  Two pairme, pairall;

  // internal methods

  void grow_bond();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
