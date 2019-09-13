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

#ifdef FIX_CLASS

FixStyle(wall/Maxwell,FixWallMaxwell)

#else

#ifndef LMP_FIX_WALL_MAXWELL_H
#define LMP_FIX_WALL_MAXWELL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallMaxwell : public Fix {
 public:
  FixWallMaxwell(class LAMMPS *, int, char **);
  virtual ~FixWallMaxwell();
  int setmask();
  void init();
  void post_integrate();

 private:
  double coeff1[6],coeff2[6];
  int seedfix;

 protected:
  int nwall;
  int wallwhich[6],wallstyle[6];
  double walltemp[6],wallvel[6][3],wallaccom[6];
  double coord0[6];
  char *varstr[6];
  int varindex[6];
  int varflag;
  double xscale,yscale,zscale;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall/gas command

Self-explanatory.

E: Cannot use fix wall/gas in periodic dimension

Self-explanatory.

E: Cannot use fix wall/reflect zlo/zhi for a 2d simulation

Self-explanatory.

E: Variable name for fix wall/gas does not exist

Self-explanatory.

E: Variable for fix wall/gas is invalid style

Only equal-style variables can be used.

W: Should not allow rigid bodies to bounce off relecting walls

LAMMPS allows this, but their dynamics are not computed correctly.

*/
