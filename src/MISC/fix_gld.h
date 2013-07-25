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

FixStyle(gld,FixGLD)

#else

#ifndef LMP_FIX_GLD_H
#define LMP_FIX_GLD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGLD : public Fix {
 public:
  FixGLD(class LAMMPS *, int, char **);
  virtual ~FixGLD();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  void reset_target(double);
  virtual void reset_dt();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void init_s_gld();
  
 protected:
  double dtv,dtf;
  double *step_respa;
  int mass_require;
  int freezeflag, zeroflag;
  double t_start, t_stop, t_target;

  int prony_terms;
  int series_type;
  double *prony_c;
  double *prony_tau;

  double **s_gld;

  class RanMars *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
