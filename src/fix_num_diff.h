/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(numdiff,FixNumDiff)

#else

#ifndef LMP_FIX_NUMDIFF_H
#define LMP_FIX_NUMDIFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNumDiff : public Fix {
 public:
  FixNumDiff(class LAMMPS *, int, char **);
  ~FixNumDiff();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_array(int, int);
  double memory_usage();

protected:
  int eflag;            // flags for energy/virial computation
  int external_force_clear;   // clear forces locally or externally

  int triclinic;              // 0 if domain is orthog, 1 if triclinic
  int pairflag;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

  double **local_forces;            // local forces from numerical difference (this might be usefull for debugging)
  double **global_forces;           // global forces from numerical difference

  void update_force(int vflag);
  void force_clear();
  // virtual void openfile(const char* filename);

 private:
  void create_groupmap();
  void displace_atom(int local_idx, int direction, int magnitude);
  void nd_force_clear(double **forces);
  void calculate_forces(int vflag);

  int ilevel_respa;
  double del;

  int igroup,groupbit;
  bigint gcount;             // number of atoms in group
  bigint flen;             // rank of dynamical matrix
  int scaleflag;
  int me;
  bigint *groupmap;
  double **temp_f;
  double *e;

  int compressed;            // 1 if dump file is written compressed, 0 no
  int binaryflag;            // 1 if dump file is written binary, 0 no
  int file_opened;           // 1 if openfile method has been called, 0 no
  int file_flag;             // 1 custom file name, 0 dynmat.dat

  FILE *fp;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix numdiff does not exist

Self-explanatory.

E: Variable name for fix numdiff does not exist

Self-explanatory.

E: Variable for fix numdiff is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix numdiff

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix numdiff

Must define an energy variable when applying a dynamic
force during minimization.

*/
