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

#include "string.h"
#include "stdlib.h"
#include "compute_reduce_region.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SUM,MINN,MAXX};
enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduceRegion::ComputeReduceRegion(LAMMPS *lmp, int narg, char **arg) :
  ComputeReduce(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

double ComputeReduceRegion::compute_one(int m)
{
  int i;

  Region *region = domain->regions[iregion];

  // invoke the appropriate attribute,compute,fix,variable
  // compute scalar quantity by summing over atom scalars
  // only include atoms in group

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int n = value2index[m];
  int j = argindex[m];

  double one;
  if (mode == SUM) one = 0.0;
  else if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;

  if (which[m] == X) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	combine(one,x[i][j]);
  } else if (which[m] == V) {
    double **v = atom->v;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	combine(one,v[i][j]);
  } else if (which[m] == F) {
    double **f = atom->f;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	combine(one,f[i][j]);
    
  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[n];
    if (!(compute->invoked_flag & INVOKED_PERATOM)) {
      compute->compute_peratom();
      compute->invoked_flag |= INVOKED_PERATOM;
    }

    if (j == 0) {
      double *compute_scalar = compute->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	  combine(one,compute_scalar[i]);
    } else {
      double **compute_vector = compute->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	  combine(one,compute_vector[i][jm1]);
    }

  // access fix fields, check if frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[n]->peratom_freq)
      error->all("Fix used in compute reduce not computed at compatible time");

    if (j == 0) {
      double *fix_scalar = modify->fix[n]->scalar_atom;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	  combine(one,fix_scalar[i]);
    } else {
      double **fix_vector = modify->fix[n]->vector_atom;
      int jm1 = j - 1;
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	  combine(one,fix_vector[i][jm1]);
    }
    
  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (nlocal > maxatom) {
      maxatom = atom->nmax;
      memory->sfree(varatom);
      varatom =	(double *) 
	memory->smalloc(maxatom*sizeof(double),"compute/reduce:varatom");
    }

    input->variable->compute_atom(n,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
	combine(one,varatom[i]);
  }

  return one;
}
