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

#ifdef ANGLE_CLASS

AngleStyle(hybrid,AngleHybrid)

#else

#ifndef LMP_ANGLE_HYBRID_H
#define LMP_ANGLE_HYBRID_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleHybrid : public Angle {
 public:
  int nstyles;                  // # of different angle styles
  Angle **styles;               // class list for each Angle style
  char **keywords;              // keyword for each Angle style

  AngleHybrid(class LAMMPS *);
  ~AngleHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);
  double memory_usage();

 private:
  int *map;                     // which style each angle type points to

  int *nanglelist;              // # of angles in sub-style anglelists
  int *maxangle;                // max # of angles sub-style lists can store
  int ***anglelist;             // anglelist for each sub-style

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

E: Angle style hybrid cannot use same angle style twice

Self-explanatory.

E: Angle style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Angle style hybrid cannot have none as an argument

Self-explanatory.

E: BondAngle coeff for hybrid angle has invalid format

No "ba" field should appear in data file entry.

E: BondBond coeff for hybrid angle has invalid format

No "bb" field should appear in data file entry.

E: Angle coeff for hybrid has invalid style

Angle style hybrid uses another angle style as one of its
coefficients.  The angle style used in the angle_coeff command or read
from a restart file is not recognized.

E: Invoked angle equil angle on angle style none

Self-explanatory.

E: Invoked angle single on angle style none

Self-explanatory.

*/
