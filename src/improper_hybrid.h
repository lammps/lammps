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

#ifdef IMPROPER_CLASS

ImproperStyle(hybrid,ImproperHybrid)

#else

#ifndef LMP_IMPROPER_HYBRID_H
#define LMP_IMPROPER_HYBRID_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperHybrid : public Improper {
 public:
  int nstyles;                  // # of different improper styles
  Improper **styles;            // class list for each Improper style
  char **keywords;              // keyword for each improper style

  ImproperHybrid(class LAMMPS *);
  ~ImproperHybrid();
  void init_style();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 private:
  int *map;                     // which style each improper type points to

  int *nimproperlist;           // # of impropers in sub-style improperlists
  int *maximproper;             // max # of impropers sub-style lists can store
  int ***improperlist;          // improperlist for each sub-style

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

E: Improper style hybrid cannot use same improper style twice

Self-explanatory.

E: Improper style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Improper style hybrid cannot have none as an argument

Self-explanatory.

E: Improper coeff for hybrid has invalid style

Improper style hybrid uses another improper style as one of its
coefficients.  The improper style used in the improper_coeff command
or read from a restart file is not recognized.

*/
