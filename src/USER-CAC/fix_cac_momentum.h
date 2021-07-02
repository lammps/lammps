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

FixStyle(cac/momentum,FixCACMomentum)

#else

#ifndef LMP_FIX_CAC_MOMENTUM_H
#define LMP_FIX_CAC_MOMENTUM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCACMomentum : public Fix {
 public:
  FixCACMomentum(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();

 protected:
  int linear,angular,rescale;
  int xflag,yflag,zflag;
  int quadrature_node_count;
  double masstotal;
  double *quadrature_weights;
  double *quadrature_abcissae;

  double shape_function(double, double, double,int,int);
  void quadrature_init();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix momentum group has no atoms

Self-explanatory.

*/
