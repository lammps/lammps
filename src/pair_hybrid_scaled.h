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
PairStyle(hybrid/scaled,PairHybridScaled);
// clang-format on
#else

#ifndef LMP_PAIR_HYBRID_SCALED_H
#define LMP_PAIR_HYBRID_SCALED_H

#include "pair_hybrid.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class PairHybridScaled : public PairHybrid {
 public:
  PairHybridScaled(class LAMMPS *);
  virtual ~PairHybridScaled();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

  void init_svector();
  void copy_svector(int, int);

 protected:
  double **fsum, **tsum;
  double *scaleval;
  int *scaleidx;
  std::vector<std::string> scalevars;
  int nmaxfsum;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

*/
