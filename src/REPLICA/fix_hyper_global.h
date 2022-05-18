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
  ~FixHyperGlobal() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup_pre_neighbor() override;
  void setup_pre_reverse(int, int) override;
  void pre_neighbor() override;
  void pre_reverse(int, int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double query(int) override;

  double memory_usage() override;

  // extra methods visible to callers

  void init_hyper() override;
  void build_bond_list(int) override;

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
