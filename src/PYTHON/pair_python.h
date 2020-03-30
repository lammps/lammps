/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.

   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(python,PairPython)

#else

#ifndef LMP_PAIR_PYTHON_H
#define LMP_PAIR_PYTHON_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPython : public Pair {
 public:
  PairPython(class LAMMPS *);
  virtual ~PairPython();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global;
  void * py_potential;
  int  * skip_types;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Could not find 'compute_force' method'

UNDOCUMENTED

E: Python 'compute_force' is not callable

UNDOCUMENTED

E: Could not find 'compute_energy' method'

UNDOCUMENTED

E: Python 'compute_energy' is not callable

UNDOCUMENTED

E: Could not create tuple for 'compute' function arguments

UNDOCUMENTED

E: Calling 'compute_force' function failed

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Python pair style requires fully qualified class name

UNDOCUMENTED

E: Loading python pair style module failure

UNDOCUMENTED

E: Could not find pair style class in module'

UNDOCUMENTED

E: Could not instantiate instance of pair style class'

UNDOCUMENTED

E: Could not find 'check_units' method'

UNDOCUMENTED

E: Python 'check_units' is not callable

UNDOCUMENTED

E: Could not create tuple for 'check_units' function arguments

UNDOCUMENTED

E: Calling 'check_units' function failed

UNDOCUMENTED

E: Could not find 'map_coeff' method'

UNDOCUMENTED

E: Python 'map_coeff' is not callable

UNDOCUMENTED

E: Could not create tuple for 'map_coeff' function arguments

UNDOCUMENTED

E: Calling 'map_coeff' function failed

UNDOCUMENTED

E: Calling 'compute_energy' function failed

UNDOCUMENTED

U: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
