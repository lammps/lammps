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
PairStyle(atm,PairATM);
// clang-format on
#else

#ifndef LMP_PAIR_ATM_H
#define LMP_PAIR_ATM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairATM : public Pair {
 public:
  PairATM(class LAMMPS *);
  virtual ~PairATM();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double cut_global, cut_triple;
  double ***nu;

  void allocate();
  void interaction_ddd(double, double, double, double, double, double *, double *, double *,
                       double *, double *, int, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal pair_style command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style ATM requires newton pair on

See the newton command.  This is a restriction to use the ATM
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
