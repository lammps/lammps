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

FixStyle(tmd,FixTMD)

#else

#ifndef LMP_FIX_TMD_H
#define LMP_FIX_TMD_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixTMD : public Fix {
 public:
  FixTMD(class LAMMPS *, int, char **);
  ~FixTMD();
  int setmask();
  void init();
  void initial_integrate(int);
  void initial_integrate_respa(int, int, int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void reset_dt();

 private:
  int me;
  int nfileevery,compressed;
  bigint previous_stat;
  FILE *fp;
  double rho_start,rho_stop,rho_old,masstotal;
  double dtv,dtf;
  double *step_respa;
  double work_lambda,work_analytical;
  double **xf,**xold;

  void readfile(char *);
  void open(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix TMD unless atom map exists

Using this fix requires the ability to lookup an atom index, which is
provided by an atom map.  An atom map does not exist (by default) for
non-molecular problems.  Using the atom_modify map command will force
an atom map to be created.

E: Cannot open fix tmd file %s

The output file for the fix tmd command cannot be opened.  Check that
the path and name are correct.

E: Fix tmd must come after integration fixes

Any fix tmd command must appear in the input script after all time
integration fixes (nve, nvt, npt).  See the fix tmd documentation for
details.

E: Incorrect format in TMD target file

Format of file read by fix tmd command is incorrect.

E: TMD target file did not list all group atoms

The target file for the fix tmd command did not list all atoms in the
fix group.

E: Cannot open gzipped file

LAMMPS is attempting to open a gzipped version of the specified file
but was unsuccessful.  Check that the path and name are correct.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
