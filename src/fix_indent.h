/* ----------------------------------------------------------------------
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

FixStyle(indent,FixIndent)

#else

#ifndef LMP_FIX_INDENT_H
#define LMP_FIX_INDENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIndent : public Fix {
 public:
  FixIndent(class LAMMPS *, int, char **);
  ~FixIndent();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  virtual void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  int istyle,scaleflag,thermo_flag,eflag_enable,side;
  double k,k3;
  char *xstr,*ystr,*zstr,*rstr,*pstr;
  int xvar,yvar,zvar,rvar,pvar;
  double xvalue,yvalue,zvalue,rvalue,pvalue;
  int indenter_flag,planeside;
  double indenter[4],indenter_all[4];
  int cdim;
  int nlevels_respa;

  void options(int, char **);
};

}

#endif
#endif
