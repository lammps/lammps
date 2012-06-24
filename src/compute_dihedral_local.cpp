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
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000
#define SMALL 0.001

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::ComputeDihedralLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute dihedral/local command");

  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,
               "Compute dihedral/local used when dihedrals are not allowed");

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
    else error->all(FLERR,"Invalid keyword in compute dihedral/local command");
  }

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::~ComputeDihedralLocal()
{
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeDihedralLocal::init()
{
  if (force->dihedral == NULL)
    error->all(FLERR,"No dihedral style is defined for compute dihedral/local");

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
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double s,c;
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
      if (pflag >= 0 && array) pbuf = &array[0][pflag];
      else pbuf = NULL;
    }
  }

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

        // phi calculation from dihedral style harmonic

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

          ax = vb1y*vb2zm - vb1z*vb2ym;
          ay = vb1z*vb2xm - vb1x*vb2zm;
          az = vb1x*vb2ym - vb1y*vb2xm;
          bx = vb3y*vb2zm - vb3z*vb2ym;
          by = vb3z*vb2xm - vb3x*vb2zm;
          bz = vb3x*vb2ym - vb3y*vb2xm;

          rasq = ax*ax + ay*ay + az*az;
          rbsq = bx*bx + by*by + bz*bz;
          rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
          rg = sqrt(rgsq);

          rginv = ra2inv = rb2inv = 0.0;
          if (rg > 0) rginv = 1.0/rg;
          if (rasq > 0) ra2inv = 1.0/rasq;
          if (rbsq > 0) rb2inv = 1.0/rbsq;
          rabinv = sqrt(ra2inv*rb2inv);

          c = (ax*bx + ay*by + az*bz)*rabinv;
          s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
          pbuf[n] = 180.0*atan2(s,c)/MY_PI;
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

double ComputeDihedralLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
