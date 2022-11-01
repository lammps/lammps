/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  friend class AtomVecDielectric;
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
  ~PairHybrid() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void setup() override;
  void finish() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;

  void modify_params(int narg, char **arg) override;
  double memory_usage() override;

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;
  void *extract(const char *, int &) override;
  void reset_dt() override;

  int check_ijtype(int, int, char *);

  void add_tally_callback(class Compute *) override;
  void del_tally_callback(class Compute *) override;
  double atom2cut(int) override;
  double radii2cut(double, double) override;

 protected:
  int nstyles;             // # of sub-styles
  Pair **styles;           // list of Pair style classes
  double *cutmax_style;    // max cutoff for each style
  char **keywords;         // style name of each Pair style
  int *multiple;           // 0 if style used once, else Mth instance

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
