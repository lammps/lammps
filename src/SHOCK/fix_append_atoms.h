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

FixStyle(append_atoms,FixAppendAtoms)

#else

#ifndef FIX_APPEND_ATOMS_H
#define FIX_APPEND_ATOMS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAppendAtoms : public Fix {
 public:
  FixAppendAtoms(class LAMMPS *, int, char **);
  ~FixAppendAtoms();
  int setmask();
  void setup(int);
  void pre_exchange();
  void initial_integrate(int);
  void post_force(int);

 private:
  int get_spatial();
  int spatflag, xloflag, xhiflag, yloflag, yhiflag, zloflag, zhiflag;
  int ranflag, tempflag, xseed, tseed;
  double ranx, rany, ranz, t_target, t_period, t_extent;
  class RanMars *randomx;
  class RanMars *randomt;
  int scaleflag, freq;
  int *basistype, nbasis;
  int advance, advance_sum;
  double size, spatlead;
  char *spatialid;
  double tfactor;
  double *gfactor1,*gfactor2;
};

}

#endif
#endif
