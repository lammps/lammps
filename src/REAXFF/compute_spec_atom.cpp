// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Labo0ratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_spec_atom.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "pair_reaxff.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{KEYWORD,COMPUTE,FIX,VARIABLE};

/* ---------------------------------------------------------------------- */

ComputeSpecAtom::ComputeSpecAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute spec/atom command");

  peratom_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // get reference to ReaxFF pair style
  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..",0));

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

    // from pair_reaxff
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
    } else if (strcmp(arg[iarg],"abo13") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo13;
    } else if (strcmp(arg[iarg],"abo14") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo14;
    } else if (strcmp(arg[iarg],"abo15") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo15;
    } else if (strcmp(arg[iarg],"abo16") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo16;
    } else if (strcmp(arg[iarg],"abo17") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo17;
    } else if (strcmp(arg[iarg],"abo18") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo18;
    } else if (strcmp(arg[iarg],"abo19") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo19;
    } else if (strcmp(arg[iarg],"abo20") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo20;
    } else if (strcmp(arg[iarg],"abo21") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo21;
    } else if (strcmp(arg[iarg],"abo22") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo22;
    } else if (strcmp(arg[iarg],"abo23") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo23;
    } else if (strcmp(arg[iarg],"abo24") == 0) {
      pack_choice[i] = &ComputeSpecAtom::pack_abo24;

    } else error->all(FLERR,"Invalid keyword in compute spec/atom command");
  }

  nmax = 0;
  vector = nullptr;
  array = nullptr;

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

  if (atom->nmax > nmax) {
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
  double bytes = (double)nmax*nvalues * sizeof(double);
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

void ComputeSpecAtom::pack_abo01(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][0];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][1];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][2];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][3];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][4];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][5];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][6];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][7];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][8];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][9];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][10];
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
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][11];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo13(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][12];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo14(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][13];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo15(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][14];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo16(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][15];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo17(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][16];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo18(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][17];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo19(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][18];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo20(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][19];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo21(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][20];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo22(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][21];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo23(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][22];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpecAtom::pack_abo24(int n)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) buf[n] = reaxff->tmpbo[i][23];
    else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */
