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

#include "math.h"
#include "string.h"
#include "compute_improper_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "improper.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000

#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ComputeImproperLocal::ComputeImproperLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute improper/local command");

  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,"Compute improper/local used when impropers are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  cflag = -1;
  nvalues = 0;

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;
    if (strcmp(arg[iarg],"chi") == 0) cflag = nvalues++;
    else error->all(FLERR,"Invalid keyword in compute improper/local command");
  }

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeImproperLocal::~ComputeImproperLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeImproperLocal::init()
{
  if (force->improper == NULL) 
    error->all(FLERR,"No improper style is defined for compute improper/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_impropers(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeImproperLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute improper info

  ncount = compute_impropers(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_impropers(1);
}

/* ----------------------------------------------------------------------
   count impropers on this proc
   only count if 2nd atom is the one storing the improper
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if flag is set, compute requested info about improper
------------------------------------------------------------------------- */

int ComputeImproperLocal::compute_impropers(int flag)
{
  int i,m,n,atom1,atom2,atom3,atom4;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c;
  double *cbuf;

  double **x = atom->x;
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (flag) {
    if (nvalues == 1) {
      if (cflag >= 0) cbuf = vector;
    } else {
      if (cflag >= 0) cbuf = &array[0][cflag];
    }
  }

  m = n = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_improper[atom2]; i++) {
      if (tag[atom2] != improper_atom2[atom2][i]) continue;
      atom1 = atom->map(improper_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(improper_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      atom4 = atom->map(improper_atom4[atom2][i]);
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      if (flag) {

	// chi calculation from improper style harmonic

	if (cflag >= 0) {
	  vb1x = x[atom1][0] - x[atom2][0];
	  vb1y = x[atom1][1] - x[atom2][1];
	  vb1z = x[atom1][2] - x[atom2][2];
	  domain->minimum_image(vb1x,vb1y,vb1z);

	  vb2x = x[atom3][0] - x[atom2][0];
	  vb2y = x[atom3][1] - x[atom2][1];
	  vb2z = x[atom3][2] - x[atom2][2];
	  domain->minimum_image(vb2x,vb2y,vb2z);

	  vb3x = x[atom4][0] - x[atom3][0];
	  vb3y = x[atom4][1] - x[atom3][1];
	  vb3z = x[atom4][2] - x[atom3][2];
	  domain->minimum_image(vb3x,vb3y,vb3z);

	  ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
	  ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
	  ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
        
	  r1 = sqrt(ss1);
	  r2 = sqrt(ss2);
	  r3 = sqrt(ss3);
        
	  c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
	  c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
	  c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

	  s1 = 1.0 - c1*c1;
	  if (s1 < SMALL) s1 = SMALL;
	  s1 = 1.0 / s1;

	  s2 = 1.0 - c2*c2;
	  if (s2 < SMALL) s2 = SMALL;
	  s2 = 1.0 / s2;

	  s12 = sqrt(s1*s2);
	  c = (c1*c2 + c0) * s12;

	  if (c > 1.0) c = 1.0;
	  if (c < -1.0) c = -1.0;
	  cbuf[n] = 180.0*acos(c)/MY_PI;
	}
	n += nvalues;
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeImproperLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/local:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeImproperLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
