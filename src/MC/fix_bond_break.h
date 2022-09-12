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
FixStyle(bond/break,FixBondBreak);
// clang-format on
#else

#ifndef LMP_FIX_BOND_BREAK_H
#define LMP_FIX_BOND_BREAK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondBreak : public Fix {
  friend class FixSRPREACT;

 public:
  FixBondBreak(class LAMMPS *, int, char **);
  ~FixBondBreak() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  void post_integrate_respa(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  int me, nprocs;
  int btype, seed;
  double cutoff, cutsq, fraction;
  int angleflag, dihedralflag, improperflag;
  bigint lastcheck;

  int breakcount, breakcounttotal;
  int nmax;
  tagint *partner, *finalpartner;
  double *distsq, *probability;

  int nbreak, maxbreak;
  tagint **broken;

  tagint *copy;

  class RanMars *random;
  int nlevels_respa;

  int commflag;
  int nbroken;
  int nangles, ndihedrals, nimpropers;

  void check_ghosts();
  void update_topology();
  void break_angles(int, tagint, tagint);
  void break_dihedrals(int, tagint, tagint);
  void break_impropers(int, tagint, tagint);
  void rebuild_special_one(int);
  int dedup(int, int, tagint *);
};

}    // namespace LAMMPS_NS

#endif
#endif
