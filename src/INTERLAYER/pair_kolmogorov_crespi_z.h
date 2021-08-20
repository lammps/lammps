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
PairStyle(kolmogorov/crespi/z,PairKolmogorovCrespiZ);
// clang-format on
#else

#ifndef LMP_PAIR_KolmogorovCrespiZ_H
#define LMP_PAIR_KolmogorovCrespiZ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairKolmogorovCrespiZ : public Pair {
 public:
  PairKolmogorovCrespiZ(class LAMMPS *);
  virtual ~PairKolmogorovCrespiZ();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);

  static constexpr int NPARAMS_PER_LINE = 11;

 protected:
  struct Param {
    double z0, C0, C2, C4, C, delta, lambda, A, S;
    double delta2inv, z06;
    int ielement, jelement;
  };
  Param *params;    // parameter set for I-J interactions

  double cut_global;
  double **offset;
  void read_file(char *);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
