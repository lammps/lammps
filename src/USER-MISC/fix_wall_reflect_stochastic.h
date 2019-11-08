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

FixStyle(wall/reflect/stochastic,FixWallReflectStochastic)

#else

#ifndef LMP_FIX_WALL_REFLECT_STOCHASTIC_H
#define LMP_FIX_WALL_REFLECT_STOCHASTIC_H


#include "random_mars.h"
#include "fix_wall_reflect.h"

namespace LAMMPS_NS {

class FixWallReflectStochastic : public FixWallReflect {
 public:
  FixWallReflectStochastic(class LAMMPS *, int, char **);
  void wall_particle(int m,int which, double coord);
  virtual ~FixWallReflectStochastic();





 protected:
    class RanMars *random;
    int seedfix;
    double walltemp[6],wallvel[6][3],wallaccom[6][3];
    int reflectionstyle;

  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall/stochastic command

Self-explanatory.

E: Cannot use fix wall/stochastic in periodic dimension

Self-explanatory.

E: Cannot use fix wall/stochastic zlo/zhi for a 2d simulation

Self-explanatory.

E: Variable name for fix wall/stochastic does not exist

Self-explanatory.

E: Variable for fix wall/stochastic is invalid style

Only equal-style variables can be used.

W: Should not allow rigid bodies to bounce off relecting walls

LAMMPS allows this, but their dynamics are not computed correctly.

*/
