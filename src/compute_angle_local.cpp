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
#include "compute_angle_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "angle.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000

enum{THETA,ENG,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeAngleLocal::ComputeAngleLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  bstyle(NULL), vvar(NULL), tstr(NULL), vstr(NULL), vlocal(NULL), alocal(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute angle/local command");

  if (atom->avec->angles_allow == 0)
    error->all(FLERR,"Compute angle/local used when angles are not allowed");

  local_flag = 1;

  // style args

  nvalues = narg - 3;
  bstyle = new int[nvalues];
  vstr = new char*[nvalues];
  vvar = new int[nvalues];

  nvalues = 0;
  tflag = 0;
  nvar = 0;

  int iarg;
  for (iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"theta") == 0) {
      bstyle[nvalues++] = THETA;
      tflag = 1;
    } else if (strcmp(arg[iarg],"eng") == 0) {
      bstyle[nvalues++] = ENG;
    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      bstyle[nvalues++] = VARIABLE;
      int n = strlen(arg[iarg]);
      vstr[nvar] = new char[n];
      strcpy(vstr[nvar],&arg[iarg][2]);
      nvar++;
    } else break;
  }

  // optional args

  setflag = 0;
  tstr = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"set") == 0) {
      setflag = 1;
      if (iarg+3 > narg) error->all(FLERR,"Illegal compute angle/local command");
      if (strcmp(arg[iarg+1],"theta") == 0) {
        delete [] tstr;
        int n = strlen(arg[iarg+2]) + 1;
        tstr = new char[n];
        strcpy(tstr,arg[iarg+2]);
        tflag = 1;
      } else error->all(FLERR,"Illegal compute angle/local command");
      iarg += 3;
    } else error->all(FLERR,"Illegal compute angle/local command");
  }

  // error check

  if (nvar) {
    if (!setflag)
      error->all(FLERR,"Compute angle/local variable requires a set variable");
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
        error->all(FLERR,"Variable name for copute angle/local does not exist");
      if (!input->variable->equalstyle(vvar[i]))
        error->all(FLERR,"Variable for compute angle/local is invalid style");
    }

    if (tstr) {
      tvar = input->variable->find(tstr);
      if (tvar < 0)
        error->all(FLERR,"Variable name for compute angle/local does not exist");
      if (!input->variable->internalstyle(tvar))
        error->all(FLERR,"Variable for compute angle/local is invalid style");
    }
  } else if (setflag)
    error->all(FLERR,"Compute angle/local set with no variable");

  // initialize output

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  nmax = 0;
  vlocal = NULL;
  alocal = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeAngleLocal::~ComputeAngleLocal()
{
  delete [] bstyle;
  for (int i = 0; i < nvar; i++) delete [] vstr[i];
  delete [] vstr;
  delete [] vvar;

  delete [] tstr;

  memory->destroy(vlocal);
  memory->destroy(alocal);
}

/* ---------------------------------------------------------------------- */

void ComputeAngleLocal::init()
{
  if (force->angle == NULL)
    error->all(FLERR,"No angle style is defined for compute angle/local");

  if (nvar) {
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
        error->all(FLERR,"Variable name for compute angle/local does not exist");
    }

    if (tstr) {
      tvar = input->variable->find(tstr);
      if (tvar < 0)
        error->all(FLERR,"Variable name for compute angle/local does not exist");
    }
  }

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_angles(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeAngleLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute angle info

  ncount = compute_angles(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_angles(1);
}

/* ----------------------------------------------------------------------
   count angles and compute angle info on this proc
   only count if 2nd atom is the one storing the angle
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if angle is deleted (type = 0), do not count
   if angle is turned off (type < 0), still count
   if flag is set, compute requested info about angle
   if angle is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeAngleLocal::compute_angles(int flag)
{
  int i,m,n,na,atom1,atom2,atom3,imol,iatom,atype,ivar;
  tagint tagprev;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,theta;
  double *ptr;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int molecular = atom->molecular;

  // loop over all atoms and their angles

  Angle *angle = force->angle;

  m = n = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == 1) na = num_angle[atom2];
    else {
      if (molindex[atom2] < 0) continue;
      imol = molindex[atom2];
      iatom = molatom[atom2];
      na = onemols[imol]->num_angle[iatom];
    }

    for (i = 0; i < na; i++) {
      if (molecular == 1) {
        if (tag[atom2] != angle_atom2[atom2][i]) continue;
        atype = angle_type[atom2][i];
        atom1 = atom->map(angle_atom1[atom2][i]);
        atom3 = atom->map(angle_atom3[atom2][i]);
      } else {
        if (tag[atom2] != onemols[imol]->angle_atom2[atom2][i]) continue;
        atype = onemols[imol]->angle_type[atom2][i];
        tagprev = tag[atom2] - iatom - 1;
        atom1 = atom->map(onemols[imol]->angle_atom1[atom2][i]+tagprev);
        atom3 = atom->map(onemols[imol]->angle_atom3[atom2][i]+tagprev);
      }

      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      if (atype == 0) continue;

      if (!flag) {
        m++;
        continue;
      }

      // theta needed by one or more outputs

      if (tflag) {
        delx1 = x[atom1][0] - x[atom2][0];
        dely1 = x[atom1][1] - x[atom2][1];
        delz1 = x[atom1][2] - x[atom2][2];
        domain->minimum_image(delx1,dely1,delz1);

        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
        r1 = sqrt(rsq1);
        
        delx2 = x[atom3][0] - x[atom2][0];
        dely2 = x[atom3][1] - x[atom2][1];
        delz2 = x[atom3][2] - x[atom2][2];
        domain->minimum_image(delx2,dely2,delz2);

        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
        r2 = sqrt(rsq2);

        // c = cosine of angle
        // theta = angle in radians

        c = delx1*delx2 + dely1*dely2 + delz1*delz2;
        c /= r1*r2;
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        theta = acos(c);
      }

      if (nvalues == 1) ptr = &vlocal[m];
      else ptr = alocal[m];

      if (nvar) {
        ivar = 0;
        if (tstr) input->variable->internal_set(tvar,theta);
      }

      for (n = 0; n < nvalues; n++) {
        switch (bstyle[n]) {
        case THETA:
          ptr[n] = 180.0*theta/MY_PI;
          break;
        case ENG:
          if (atype > 0) ptr[n] = angle->single(atype,atom1,atom2,atom3);
          else ptr[n] = 0.0;
          break;
        case VARIABLE:
          ptr[n] = input->variable->compute_equal(vvar[ivar]);
          ivar++;
          break;
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAngleLocal::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"angle/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"angle/local:array_local");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeAngleLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
