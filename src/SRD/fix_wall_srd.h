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

FixStyle(wall/srd,FixWallSRD)

#else

#ifndef LMP_FIX_WALL_SRD_H
#define LMP_FIX_WALL_SRD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallSRD : public Fix {
 public:
  int nwall,varflag,overlap;
  int wallwhich[6];
  double xwall[6],xwallhold[6],vwall[6];
  double **fwall;

  FixWallSRD(class LAMMPS *, int, char **);
  ~FixWallSRD();
  int setmask();
  void init();
  double compute_array(int, int);

  void wall_params(int);

 private:
  int wallstyle[6];
  double coord0[6];
  char *varstr[6];
  int varindex[6];

  double dt;
  double xwalllast[6];
  int laststep;
  
  double **fwall_all;
  int force_flag;
};

}

#endif
#endif
