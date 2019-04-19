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

#ifdef BOND_CLASS

BondStyle(hybrid,BondHybrid)

#else

#ifndef LMP_BOND_HYBRID_H
#define LMP_BOND_HYBRID_H

#include <cstdio>
#include "bond.h"

namespace LAMMPS_NS {

class BondHybrid : public Bond {
  friend class Force;

 public:
  int nstyles;                  // # of different bond styles
  Bond **styles;                // class list for each Bond style
  char **keywords;              // keyword for each Bond style

  BondHybrid(class LAMMPS *);
  ~BondHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int, double &);
  double memory_usage();

 private:
  int *map;                     // which style each bond type points to
  int has_quartic;              // which style, if any is a quartic bond style
  int *nbondlist;               // # of bonds in sub-style bondlists
  int *maxbond;                 // max # of bonds sub-style lists can store
  int ***bondlist;              // bondlist for each sub-style

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Bond style hybrid cannot use same bond style twice

Self-explanatory.

E: Bond style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Bond style hybrid cannot have none as an argument

Self-explanatory.

E: Bond coeff for hybrid has invalid style

Bond style hybrid uses another bond style as one of its coefficients.
The bond style used in the bond_coeff command or read from a restart
file is not recognized.

E: Invoked bond equil distance on bond style none

Self-explanatory.

E: Invoked bond single on bond style none

Self-explanatory.

*/
