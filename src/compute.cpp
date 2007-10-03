/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory. // 
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "lattice.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Compute::Compute(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  if (narg < 3) error->all("Illegal compute command");

  // compute ID, group, and style

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  vector = NULL;
  scalar_atom = NULL;
  vector_atom = NULL;
  
  scalar_flag = vector_flag = peratom_flag = 0;
  tempflag = pressflag = 0;
  npre = 0;
  id_pre = NULL;
  comm_forward = comm_reverse = 0;

  // set modify defaults

  extra_dof = 3;
  dynamic = 0;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  delete [] id;
  delete [] style;

  for (int i = 0; i < npre; i++) delete [] id_pre[i];
  delete [] id_pre;
}

/* ---------------------------------------------------------------------- */

void Compute::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal compute_modify command");

  int iarg = 0;  while (iarg < narg) {
    if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute_modify command");
      extra_dof = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) dynamic = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) dynamic = 1;
      else error->all("Illegal compute_modify command");
      iarg += 2;
    } else error->all("Illegal compute_modify command");
  }
}
