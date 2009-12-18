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
#include "compute_dihedral_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "dihedral.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::ComputeDihedralLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute dihedral/local command");

  if (atom->avec->dihedrals_allow == 0)
    error->all("Compute dihedral/local used when dihedrals are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pflag = -1;
  nvalues = 0;

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;
    if (strcmp(arg[iarg],"phi") == 0) pflag = nvalues++;
    else error->all("Invalid keyword in compute dihedral/local command");
  }

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::~ComputeDihedralLocal()
{
  memory->sfree(vector);
  memory->destroy_2d_double_array(array);
}

/* ---------------------------------------------------------------------- */

void ComputeDihedralLocal::init()
{
  if (force->dihedral == NULL) 
    error->all("No dihedral style is defined for compute dihedral/local");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_dihedrals(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeDihedralLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute dihedral info

  ncount = compute_dihedrals(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_dihedrals(1);
}

/* ----------------------------------------------------------------------
   count dihedrals on this proc
   only count if 2nd atom is the one storing the dihedral
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if flag is set, compute requested info about dihedral
------------------------------------------------------------------------- */

int ComputeDihedralLocal::compute_dihedrals(int flag)
{
  int i,m,n,atom1,atom2,atom3,atom4;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sin2,sc1,sc2,s1,s2,s12,c;
  double *pbuf;

  double **x = atom->x;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (flag) {
    if (nvalues == 1) {
      if (pflag >= 0) pbuf = vector;
    } else {
      if (pflag >= 0) pbuf = &array[0][pflag];
    }
  }

  double PI = 4.0*atan(1.0);

  m = n = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_dihedral[atom2]; i++) {
      if (tag[atom2] != dihedral_atom2[atom2][i]) continue;
      atom1 = atom->map(dihedral_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(dihedral_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      atom4 = atom->map(dihedral_atom4[atom2][i]);
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      if (flag) {

	// phi calculation from dihedral style OPLS

	if (pflag >= 0) {
	  vb1x = x[atom1][0] - x[atom2][0];
	  vb1y = x[atom1][1] - x[atom2][1];
	  vb1z = x[atom1][2] - x[atom2][2];
	  domain->minimum_image(vb1x,vb1y,vb1z);

	  vb2x = x[atom3][0] - x[atom2][0];
	  vb2y = x[atom3][1] - x[atom2][1];
	  vb2z = x[atom3][2] - x[atom2][2];
	  domain->minimum_image(vb2x,vb2y,vb2z);

	  vb2xm = -vb2x;
	  vb2ym = -vb2y;
	  vb2zm = -vb2z;
	  domain->minimum_image(vb2xm,vb2ym,vb2zm);

	  vb3x = x[atom4][0] - x[atom3][0];
	  vb3y = x[atom4][1] - x[atom3][1];
	  vb3z = x[atom4][2] - x[atom3][2];
	  domain->minimum_image(vb3x,vb3y,vb3z);

	  sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
	  sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
	  sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

	  rb1 = sqrt(sb1);
	  rb3 = sqrt(sb3);

	  c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

	  b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
	  b1mag = sqrt(b1mag2);
	  b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
	  b2mag = sqrt(b2mag2);
	  b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
	  b3mag = sqrt(b3mag2);
	  
	  ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
	  r12c1 = 1.0 / (b1mag*b2mag);
	  c1mag = ctmp * r12c1;
	  
	  ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
	  r12c2 = 1.0 / (b2mag*b3mag);
	  c2mag = ctmp * r12c2;

	  sin2 = MAX(1.0 - c1mag*c1mag,0.0);
	  sc1 = sqrt(sin2);
	  if (sc1 < SMALL) sc1 = SMALL;
	  sc1 = 1.0/sc1;

	  sin2 = MAX(1.0 - c2mag*c2mag,0.0);
	  sc2 = sqrt(sin2);
	  if (sc2 < SMALL) sc2 = SMALL;
	  sc2 = 1.0/sc2;

	  s1 = sc1 * sc1;
	  s2 = sc2 * sc2;
	  s12 = sc1 * sc2;
	  c = (c0 + c1mag*c2mag) * s12;

	  if (c > 1.0) c = 1.0;
	  if (c < -1.0) c = -1.0;
	  pbuf[n] = 180.0*acos(c)/PI;
	}
	n += nvalues;
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeDihedralLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->sfree(vector);
    vector = (double *) memory->smalloc(nmax*sizeof(double),
					"bond/local:vector");
    vector_local = vector;
  } else {
    memory->destroy_2d_double_array(array);
    array = memory->create_2d_double_array(nmax,nvalues,
					   "bond/local:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeDihedralLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
