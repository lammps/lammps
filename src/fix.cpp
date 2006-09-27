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

#include "string.h"
#include "fix.h"
#include "group.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Fix::Fix(int narg, char **arg)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  restart_global = 0;
  restart_peratom = 0;
  force_reneighbor = 0;
  thermo_print = 0;
  thermo_energy = 0;
  nevery = 1;

  neigh_half_once = neigh_half_every = 0;
  neigh_full_once = neigh_full_every = 0;

  // mask settings - same as define settings in modify.cpp

  INITIAL_INTEGRATE = 1;
  PRE_EXCHANGE = 2;
  PRE_NEIGHBOR = 4;
  POST_FORCE = 8;
  FINAL_INTEGRATE = 16;
  END_OF_STEP = 32;
  THERMO = 64;
  INITIAL_INTEGRATE_RESPA = 128;
  POST_FORCE_RESPA = 256;
  FINAL_INTEGRATE_RESPA = 512;
  MIN_POST_FORCE = 1024;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  delete [] id;
  delete [] style;
}

/* ----------------------------------------------------------------------
   process params common to all fixes here
   if unknown param, call modify_param specific to the fix
------------------------------------------------------------------------- */

void Fix::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal fix_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"thermo") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermo_print = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) thermo_print = 1;
      else error->all("Illegal fix_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermo_energy = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) thermo_energy = 1;
      else error->all("Illegal fix_modify command");
      iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all("Illegal fix_modify command");
      iarg += n;
    }
  }
}
