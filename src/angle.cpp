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
#include "angle.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "suffix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

Angle::Angle(LAMMPS *lmp) : Pointers(lmp)
{
  energy = 0.0;
  writedata = 1;

  allocated = 0;
  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = 0;
  eatom = NULL;
  vatom = NULL;
  setflag = NULL;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;
}

/* ---------------------------------------------------------------------- */

Angle::~Angle()
{
  memory->destroy(eatom);
  memory->destroy(vatom);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Angle::init()
{
  if (!allocated && atom->nangletypes) 
    error->all(FLERR,"Angle coeffs are not set");
  for (int i = 1; i <= atom->nangletypes; i++)
    if (setflag[i] == 0) error->all(FLERR,"All angle coeffs are not set");

  init_style();
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void Angle::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,comm->nthreads*maxeatom,"bond:eatom");
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,comm->nthreads*maxvatom,6,"bond:vatom");
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

void Angle::ev_tally(int i, int j, int k, int nlocal, int newton_bond,
                     double eangle, double *f1, double *f3,
                     double delx1, double dely1, double delz1,
                     double delx2, double dely2, double delz2)
{
  double eanglethird,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy += eangle;
      else {
        eanglethird = THIRD*eangle;
        if (i < nlocal) energy += eanglethird;
        if (j < nlocal) energy += eanglethird;
        if (k < nlocal) energy += eanglethird;
      }
    }
    if (eflag_atom) {
      eanglethird = THIRD*eangle;
      if (newton_bond || i < nlocal) eatom[i] += eanglethird;
      if (newton_bond || j < nlocal) eatom[j] += eanglethird;
      if (newton_bond || k < nlocal) eatom[k] += eanglethird;
    }
  }

  if (vflag_either) {
    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += THIRD*v[0];
          virial[1] += THIRD*v[1];
          virial[2] += THIRD*v[2];
          virial[3] += THIRD*v[3];
          virial[4] += THIRD*v[4];
          virial[5] += THIRD*v[5];
        }
        if (j < nlocal) {
          virial[0] += THIRD*v[0];
          virial[1] += THIRD*v[1];
          virial[2] += THIRD*v[2];
          virial[3] += THIRD*v[3];
          virial[4] += THIRD*v[4];
          virial[5] += THIRD*v[5];
        }
        if (k < nlocal) {
          virial[0] += THIRD*v[0];
          virial[1] += THIRD*v[1];
          virial[2] += THIRD*v[2];
          virial[3] += THIRD*v[3];
          virial[4] += THIRD*v[4];
          virial[5] += THIRD*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += THIRD*v[0];
        vatom[i][1] += THIRD*v[1];
        vatom[i][2] += THIRD*v[2];
        vatom[i][3] += THIRD*v[3];
        vatom[i][4] += THIRD*v[4];
        vatom[i][5] += THIRD*v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += THIRD*v[0];
        vatom[j][1] += THIRD*v[1];
        vatom[j][2] += THIRD*v[2];
        vatom[j][3] += THIRD*v[3];
        vatom[j][4] += THIRD*v[4];
        vatom[j][5] += THIRD*v[5];
      }
      if (newton_bond || k < nlocal) {
        vatom[k][0] += THIRD*v[0];
        vatom[k][1] += THIRD*v[1];
        vatom[k][2] += THIRD*v[2];
        vatom[k][3] += THIRD*v[3];
        vatom[k][4] += THIRD*v[4];
        vatom[k][5] += THIRD*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double Angle::memory_usage()
{
  double bytes = comm->nthreads*maxeatom * sizeof(double);
  bytes += comm->nthreads*maxvatom*6 * sizeof(double);
  return bytes;
}
