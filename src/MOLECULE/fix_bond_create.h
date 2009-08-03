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

#ifndef FIX_BOND_CREATE_H
#define FIX_BOND_CREATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondCreate : public Fix {
 public:
  FixBondCreate(class LAMMPS *, int, char **);
  ~FixBondCreate();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();
  void post_integrate_respa(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me;
  int iatomtype,jatomtype;
  int btype,seed;
  int imaxbond,jmaxbond;
  int inewtype,jnewtype;
  double cutsq,fraction;

  int createcount,createcounttotal;   // bond formation stats

  int nmax;
  int *bondcount;        // count of created bonds this atom is part of
  int *partner;          // ID of preferred atom for this atom to bond to
  double *distsq;        // distance to preferred bond partner
  double *probability;   // random # to use in decision to form bond

  class RanMars *random;
  class NeighList *list;
  int countflag,commflag;
  int nlevels_respa;
};

}

#endif
