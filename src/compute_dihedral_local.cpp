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
#include "compute_dihedral_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "dihedral.h"
#include "input.h"
#include "variable.h"

#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000
#define SMALL 0.001

enum{PHI,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::ComputeDihedralLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  bstyle(NULL), vvar(NULL), pstr(NULL), vstr(NULL), vlocal(NULL), alocal(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute dihedral/local command");

  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR,
               "Compute dihedral/local used when dihedrals are not allowed");

  local_flag = 1;

  // style args

  nvalues = narg - 3;
  bstyle = new int[nvalues];
  vstr = new char*[nvalues];
  vvar = new int[nvalues];

  nvalues = 0;
  nvar = 0;

  int iarg;
  for (iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"phi") == 0) {
      bstyle[nvalues++] = PHI;
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
  pstr = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"set") == 0) {
      setflag = 1;
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal compute dihedral/local command");
      if (strcmp(arg[iarg+1],"phi") == 0) {
        delete [] pstr;
        int n = strlen(arg[iarg+2]) + 1;
        pstr = new char[n];
        strcpy(pstr,arg[iarg+2]);
      } else error->all(FLERR,"Illegal compute dihedral/local command");
      iarg += 3;
    } else error->all(FLERR,"Illegal compute dihedral/local command");
  }

  // error check

  if (nvar) {
    if (!setflag)
      error->all(FLERR,"Compute dihedral/local variable requires a set variable");
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
	error->all(FLERR,
                   "Variable name for copute dihedral/local does not exist");
      if (!input->variable->equalstyle(vvar[i]))
	error->all(FLERR,"Variable for compute dihedral/local is invalid style");
    }

    if (pstr) {
      pvar = input->variable->find(pstr);
      if (pvar < 0)
        error->all(FLERR,
                   "Variable name for compute dihedral/local does not exist");
      if (!input->variable->internalstyle(pvar))
        error->all(FLERR,"Variable for compute dihedral/local is invalid style");
    }
  } else if (setflag)
    error->all(FLERR,"Compute dihedral/local set with no variable");

  // initialize output

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  nmax = 0;
  vlocal = NULL;
  alocal = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDihedralLocal::~ComputeDihedralLocal()
{
  delete [] bstyle;
  for (int i = 0; i < nvar; i++) delete [] vstr[i];
  delete [] vstr;
  delete [] vvar;

  delete [] pstr;

  memory->destroy(vlocal);
  memory->destroy(alocal);
}

/* ---------------------------------------------------------------------- */

void ComputeDihedralLocal::init()
{
  if (force->dihedral == NULL)
    error->all(FLERR,"No dihedral style is defined for compute dihedral/local");

  if (nvar) {
    for (int i = 0; i < nvar; i++) {
      vvar[i] = input->variable->find(vstr[i]);
      if (vvar[i] < 0)
	error->all(FLERR,
                   "Variable name for compute dihedral/local does not exist");
    }

    if (pstr) {
      pvar = input->variable->find(pstr);
      if (pvar < 0)
        error->all(FLERR,
                   "Variable name for compute dihedral/local does not exist");
    }
  }

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
  int i,m,n,nd,atom1,atom2,atom3,atom4,imol,iatom,ivar;
  tagint tagprev;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,ra2inv,rb2inv,rabinv;
  double s,c,phi;
  double *ptr;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int *mask = atom->mask;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int molecular = atom->molecular;

  // loop over all atoms and their dihedrals

  m = n = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;

    if (molecular == 1) nd = num_dihedral[atom2];
    else {
      if (molindex[atom2] < 0) continue;
      imol = molindex[atom2];
      iatom = molatom[atom2];
      nd = onemols[imol]->num_dihedral[iatom];
    }

    for (i = 0; i < nd; i++) {
      if (molecular == 1) {
        if (tag[atom2] != dihedral_atom2[atom2][i]) continue;
        atom1 = atom->map(dihedral_atom1[atom2][i]);
        atom3 = atom->map(dihedral_atom3[atom2][i]);
        atom4 = atom->map(dihedral_atom4[atom2][i]);
      } else {
        if (tag[atom2] != onemols[imol]->dihedral_atom2[atom2][i]) continue;
        tagprev = tag[atom2] - iatom - 1;
        atom1 = atom->map(onemols[imol]->dihedral_atom1[atom2][i]+tagprev);
        atom3 = atom->map(onemols[imol]->dihedral_atom3[atom2][i]+tagprev);
        atom4 = atom->map(onemols[imol]->dihedral_atom4[atom2][i]+tagprev);
      }

      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      if (!flag) {
        m++;
        continue;
      }

      // phi calculation from dihedral style harmonic

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

      ra2inv = rb2inv = 0.0;
      if (rasq > 0) ra2inv = 1.0/rasq;
      if (rbsq > 0) rb2inv = 1.0/rbsq;
      rabinv = sqrt(ra2inv*rb2inv);

      c = (ax*bx + ay*by + az*bz)*rabinv;
      s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      phi = atan2(s,c);

      if (nvalues == 1) ptr = &vlocal[m];
      else ptr = alocal[m];

      if (nvar) {
	ivar = 0;
	if (pstr) input->variable->internal_set(pvar,phi);
      }

      for (n = 0; n < nvalues; n++) {
	switch (bstyle[n]) {
	case PHI:
	  ptr[n] = 180.0*phi/MY_PI;
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

void ComputeDihedralLocal::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"dihedral/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"dihedral/local:array_local");
    array_local = alocal;
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
