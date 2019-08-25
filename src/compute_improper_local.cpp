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

#include <cmath>
#include <cstring>
#include "compute_improper_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
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
  Compute(lmp, narg, arg),
  vlocal(NULL), alocal(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute improper/local command");

  if (atom->avec->impropers_allow == 0)
    error->all(FLERR,
               "Compute improper/local used when impropers are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  cflag = -1;
  nvalues = 0;

  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"chi") == 0) cflag = nvalues++;
    else error->all(FLERR,"Invalid keyword in compute improper/local command");
  }

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  nmax = 0;
  vlocal = NULL;
  alocal = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeImproperLocal::~ComputeImproperLocal()
{
  memory->destroy(vlocal);
  memory->destroy(alocal);
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
  int i,m,n,ni,atom1,atom2,atom3,atom4,imol,iatom;
  tagint tagprev;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c;
  double *cbuf;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_improper = atom->num_improper;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int molecular = atom->molecular;

  if (flag) {
    if (nvalues == 1) {
      if (cflag >= 0) cbuf = vlocal;
    } else {
      if (cflag >= 0 && alocal) cbuf = &alocal[0][cflag];
      else cbuf = NULL;
    }
  }

  m = n = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == 1) ni = num_improper[atom2];
    else {
      if (molindex[atom2] < 0) continue;
      imol = molindex[atom2];
      iatom = molatom[atom2];
      ni = onemols[imol]->num_improper[iatom];
    }

    for (i = 0; i < ni; i++) {
      if (molecular == 1) {
        if (tag[atom2] != improper_atom2[atom2][i]) continue;
        atom1 = atom->map(improper_atom1[atom2][i]);
        atom3 = atom->map(improper_atom3[atom2][i]);
        atom4 = atom->map(improper_atom4[atom2][i]);
      } else {
        if (tag[atom2] != onemols[imol]->improper_atom2[atom2][i]) continue;
        tagprev = tag[atom2] - iatom - 1;
        atom1 = atom->map(onemols[imol]->improper_atom1[atom2][i]+tagprev);
        atom3 = atom->map(onemols[imol]->improper_atom3[atom2][i]+tagprev);
        atom4 = atom->map(onemols[imol]->improper_atom4[atom2][i]+tagprev);
      }

      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
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
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"improper/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"improper/local:array_local");
    array_local = alocal;
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
