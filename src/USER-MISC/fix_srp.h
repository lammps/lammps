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

FixStyle(SRP,FixSRP)

#else

#ifndef LMP_FIX_SRP_H
#define LMP_FIX_SRP_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixSRP : public Fix {
 public:
  FixSRP(class LAMMPS *, int, char **);
  ~FixSRP();
  int setmask();
  void init();

  void pre_exchange();
  void setup_pre_force(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_border(int, int *, double *);
  int unpack_border(int, int, double *);
  void post_run(); 

  int pack_restart(int, double*);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);

  double **array;

 private:
  double xone[3];
  int btype;
  int bptype;
  int setup;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
