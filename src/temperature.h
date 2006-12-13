/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "lammps.h"

class Temperature : public LAMMPS {
 public:
  char *id,*style;
  int igroup,groupbit;
  double t_total,dof;
  double ke_tensor[6];

  Temperature(int, char **);
  virtual ~Temperature();
  void modify_params(int, char **);
  void count_fix();
  virtual void init() = 0;
  virtual double compute() = 0;
  virtual void tensor() = 0;

 protected:
  double tfactor;
  int extra_dof,fix_dof;
  int scaleflag,dynamic;
  double xscale,yscale,zscale;

  void options(int, char **);
};

#endif
