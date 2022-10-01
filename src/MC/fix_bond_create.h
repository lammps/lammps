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
FixStyle(bond/create,FixBondCreate);
// clang-format on
#else

#ifndef LMP_FIX_BOND_CREATE_H
#define LMP_FIX_BOND_CREATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondCreate : public Fix {
  friend class FixSRPREACT;

 public:
  FixBondCreate(class LAMMPS *, int, char **);
  ~FixBondCreate() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void post_integrate() override;
  void post_integrate_respa(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  double compute_vector(int) override;
  double memory_usage() override;

 protected:
  int me;
  int iatomtype, jatomtype;
  int btype, seed;
  int imaxbond, jmaxbond;
  int inewtype, jnewtype;
  int constrainflag, constrainpass;
  double amin, amax;
  double cutsq, fraction;
  int atype, dtype, itype;
  int angleflag, dihedralflag, improperflag;

  int overflow;
  tagint lastcheck;

  int *bondcount;
  int createcount, createcounttotal;
  int nmax;
  tagint *partner, *finalpartner;
  double *distsq, *probability;

  int ncreate, maxcreate;
  tagint **created;

  tagint *copy;

  class RanMars *random;
  class NeighList *list;

  int countflag, commflag;
  int nlevels_respa;
  int nangles, ndihedrals, nimpropers;

  void check_ghosts();
  void update_topology();
  void rebuild_special_one(int);
  void create_angles(int);
  void create_dihedrals(int);
  void create_impropers(int);
  int dedup(int, int, tagint *);

  virtual int constrain(int, int, double, double) { return 1; }
};

}    // namespace LAMMPS_NS

#endif
#endif
