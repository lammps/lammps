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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hybrid,PairHybrid);
// clang-format on
#else

#ifndef LMP_PAIR_HYBRID_H
#define LMP_PAIR_HYBRID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairHybrid : public Pair {
  friend class ComputeSpin;
  friend class FixGPU;
  friend class FixIntel;
  friend class FixNVESpin;
  friend class FixOMP;
  friend class Force;
  friend class Info;
  friend class Neighbor;
  friend class PairDeprecated;
  friend class Respa;
  friend class Scafacos;

 public:
  PairHybrid(class LAMMPS *);
  virtual ~PairHybrid();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void setup();
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  void modify_params(int narg, char **arg);
  double memory_usage();

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
  void *extract(const char *, int &);
  void reset_dt();

  int check_ijtype(int, int, char *);

  virtual void add_tally_callback(class Compute *);
  virtual void del_tally_callback(class Compute *);
  double atom2cut(int);
  double radii2cut(double, double);

 protected:
  int nstyles;        // # of sub-styles
  Pair **styles;      // list of Pair style classes
  char **keywords;    // style name of each Pair style
  int *multiple;      // 0 if style used once, else Mth instance

  int outerflag;    // toggle compute() when invoked by outer()
  int respaflag;    // 1 if different substyles are assigned to
                    // different r-RESPA levels

  int **nmap;               // # of sub-styles itype,jtype points to
  int ***map;               // list of sub-styles itype,jtype points to
  double **special_lj;      // list of per style LJ exclusion factors
  double **special_coul;    // list of per style Coulomb exclusion factors
  int *compute_tally;       // list of on/off flags for tally computes

  void allocate();
  void flags();

  virtual void init_svector();
  virtual void copy_svector(int, int);

  void modify_special(int, int, char **);
  double *save_special();
  void set_special(int);
  void restore_special(double *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Pair style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Pair style hybrid cannot have none as an argument

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot yet use pair hybrid with Kokkos

This feature is not yet supported.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

E: Pair hybrid sub-style is not used

No pair_coeff command used a sub-style specified in the pair_style
command.

E: Pair_modify special setting for pair hybrid incompatible with global special_bonds setting

Cannot override a setting of 0.0 or 1.0 or change a setting between
0.0 and 1.0.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Invoked pair single on pair style none

A command (e.g. a dump) attempted to invoke the single() function on a
pair style none, which is illegal.  You are probably attempting to
compute per-atom quantities with an undefined pair style.

E: Pair hybrid sub-style does not support single call

You are attempting to invoke a single() call on a pair style
that doesn't support it.

E: Pair hybrid single calls do not support per sub-style special bond values

Self-explanatory.

E: Unknown pair_modify hybrid sub-style

The choice of sub-style is unknown.

E: Coulomb cutoffs of pair hybrid sub-styles do not match

If using a Kspace solver, all Coulomb cutoffs of long pair styles must
be the same.

*/
