/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Labo0ratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "compute_spec_atom.h"
#include "math_extra.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#include "reaxc_defs.h"
#include "reaxc_types.h"
#include "pair_reax_c.h"

using namespace LAMMPS_NS;

enum{KEYWORD,COMPUTE,FIX,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeSpecAtom::ComputeSpecAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute reax/c/atom command");

  peratom_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // Initiate reaxc
  reaxc = (PairReaxC *) force->pair_match("reax/c",1);

  pack_choice = new FnPtrPack[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    // standard lammps attributes
    if (strcmp(arg[iarg],"q") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_q;

    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_x;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_y;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_z;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_vx;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_vy;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_vz;

    // from pair_reax_c
    } else if (strcmp(arg[iarg],"jid01") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid01;
    } else if (strcmp(arg[iarg],"jid02") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid02;
    } else if (strcmp(arg[iarg],"jid03") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid03;
    } else if (strcmp(arg[iarg],"jid04") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid04;
    } else if (strcmp(arg[iarg],"jid05") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid05;
    } else if (strcmp(arg[iarg],"jid06") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid06;
    } else if (strcmp(arg[iarg],"jid07") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid07;
    } else if (strcmp(arg[iarg],"jid08") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid08;
    } else if (strcmp(arg[iarg],"jid09") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid09;
    } else if (strcmp(arg[iarg],"jid10") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid10;
    } else if (strcmp(arg[iarg],"jid11") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid11;
    } else if (strcmp(arg[iarg],"jid12") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_jid12;
    } else if (strcmp(arg[iarg],"abo01") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo01;
    } else if (strcmp(arg[iarg],"abo02") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo02;
    } else if (strcmp(arg[iarg],"abo03") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo03;
    } else if (strcmp(arg[iarg],"abo04") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo04;
    } else if (strcmp(arg[iarg],"abo05") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo05;
    } else if (strcmp(arg[iarg],"abo06") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo06;
    } else if (strcmp(arg[iarg],"abo07") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo07;
    } else if (strcmp(arg[iarg],"abo08") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo08;
    } else if (strcmp(arg[iarg],"abo09") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo09;
    } else if (strcmp(arg[iarg],"abo10") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo10;
    } else if (strcmp(arg[iarg],"abo11") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo11;
    } else if (strcmp(arg[iarg],"abo12") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo12;

    } else error->all(FLERR,"Invalid keyword in compute reax/c/atom command");
  }

  nmax = 0;
  vector = NULL;
  array = NULL;

}

/* ---------------------------------------------------------------------- */

ComputeSpecAtom::~ComputeSpecAtom()
{
  delete [] pack_choice;
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow vector or array if necessary

  if (atom->nlocal > nmax) {
    nmax = atom->nmax;
    if (nvalues == 1) {
      memory->destroy(vector);
      memory->create(vector,nmax,"property/atom:vector");
      vector_atom = vector;
    } else {
      memory->destroy(array);
      memory->create(array,nmax,nvalues,"property/atom:array");
      array_atom = array;
    }
  }

  // fill vector or array with per-atom values

  if (nvalues == 1) {
    buf = vector;
    (this->*pack_choice[0])(0);
  } else {
    if (nmax > 0) {
      buf = &array[0][0];
      for (int n = 0; n < nvalues; n++)
        (this->*pack_choice[n])(n);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSpecAtom::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/atom can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_q(int n)
{
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = q[i];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_x(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = x[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_y(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = x[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_z(int n)
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = x[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_vx(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = v[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_vy(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = v[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_vz(int n)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = v[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid01(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid02(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid03(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid04(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][3];
    else buf[n] = 0.0;
    n += nvalues;
  }
}


/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid05(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][4];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid06(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][5];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid07(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][6];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid08(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][7];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid09(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][8];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid10(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][9];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid11(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][10];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo01(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][0];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo02(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][1];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo03(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][2];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo04(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][3];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo05(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][4];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo06(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][5];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo07(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][6];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo08(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][7];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo09(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][8];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo10(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][9];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo11(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][10];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_jid12(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpid[i][11];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo12(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxc->tmpbo[i][11];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

