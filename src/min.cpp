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

#include "stdlib.h"
#include "string.h"
#include "min.h"
#include "error.h"

using namespace LAMMPS_NS;

#define SCAN   0   // same as in min_cg.cpp
#define SECANT 1

/* ---------------------------------------------------------------------- */

Min::Min(LAMMPS *lmp) : Pointers(lmp)
{
  linestyle = SECANT;
  dmin = 1.0e-5;
  dmax = 0.1;
  lineiter = 10;
}

/* ---------------------------------------------------------------------- */

void Min::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal min_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linestyle") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      if (strcmp(arg[iarg+1],"scan") == 0) linestyle = SCAN;
      else if (strcmp(arg[iarg+1],"secant") == 0) linestyle = SECANT;
      else error->all("Illegal min_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dmin") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      dmin = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dmax") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      dmax = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"lineiter") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      lineiter = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal min_modify command");
  }
}
