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
/* ----------------------------------------------------------------------
   Contributing authors: Eugen Rozic (University College London)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(cosine/squared, PairCosineSquared)

#else

#ifndef LMP_PAIR_LJ_COS_SQ_H
#define LMP_PAIR_LJ_COS_SQ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCosineSquared : public Pair {
 public:
  PairCosineSquared(class LAMMPS *);
  virtual ~PairCosineSquared();
  void settings(int, char **);
  void coeff(int, char **);
  // void init_style();
  double init_one(int, int);
  void modify_params(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual void compute(int, int);
  double single(int, int, int, int, double, double, double, double &);
  // void *extract(const char *, int &);

/* RESPA stuff not implemented...
  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
*/

 protected:
  double cut_global;
  double **epsilon, **sigma, **w, **cut;
  int **wcaflag;
  double **lj12_e, **lj6_e, **lj12_f, **lj6_f;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Mixing not supported in pair_style cosine/squared

Self-explanatory. All coefficients need to be specified explicitly.

E: pair_modify mix not supported for pair_style cosine/squared

Same as above, only when calling "pair_modify" command

W: pair_modify shift/tail is meaningless for pair_style cosine/squared

This style by definition gets to zero at cutoff distance, so there is nothing
to shift and there is no tail contribution

*/

