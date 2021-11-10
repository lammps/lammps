// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "improper.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "suffix.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Improper::Improper(LAMMPS *lmp) : Pointers(lmp)
{
  energy = 0.0;
  writedata = 0;

  allocated = 0;
  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = maxcvatom = 0;
  eatom = nullptr;
  vatom = nullptr;
  cvatom = nullptr;
  setflag = nullptr;
  centroidstressflag = CENTROID_AVAIL;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Improper::~Improper()
{
  if (copymode) return;

  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(cvatom);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Improper::init()
{
  if (!allocated && atom->nimpropertypes)
    error->all(FLERR,"Improper coeffs are not set");
  for (int i = 1; i <= atom->nimpropertypes; i++)
    if (setflag[i] == 0) error->all(FLERR,"All improper coeffs are not set");

  init_style();
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for bitwise settings of eflag/vflag
   set the following flags, values are otherwise set to 0:
     evflag       != 0 if any bits of eflag or vflag are set
     eflag_global != 0 if ENERGY_GLOBAL bit of eflag set
     eflag_atom   != 0 if ENERGY_ATOM bit of eflag set
     eflag_either != 0 if eflag_global or eflag_atom is set
     vflag_global != 0 if VIRIAL_PAIR or VIRIAL_FDOTR bit of vflag set
     vflag_atom   != 0 if VIRIAL_ATOM bit of vflag set
     vflag_atom   != 0 if VIRIAL_CENTROID bit of vflag set
                       and centroidstressflag != CENTROID_AVAIL
     cvflag_atom  != 0 if VIRIAL_CENTROID bit of vflag set
                       and centroidstressflag = CENTROID_AVAIL
     vflag_either != 0 if any of vflag_global, vflag_atom, cvflag_atom is set
------------------------------------------------------------------------- */

void Improper::ev_setup(int eflag, int vflag, int alloc)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  vflag_atom = vflag & VIRIAL_ATOM;
  if (vflag & VIRIAL_CENTROID && centroidstressflag != CENTROID_AVAIL)
    vflag_atom = 1;
  cvflag_atom = 0;
  if (vflag & VIRIAL_CENTROID && centroidstressflag == CENTROID_AVAIL)
    cvflag_atom = 1;
  vflag_either = vflag_global || vflag_atom || cvflag_atom;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom,comm->nthreads*maxeatom,"improper:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom,comm->nthreads*maxvatom,6,"improper:vatom");
    }
  }
  if (cvflag_atom && atom->nmax > maxcvatom) {
    maxcvatom = atom->nmax;
    if (alloc) {
      memory->destroy(cvatom);
      memory->create(cvatom,comm->nthreads*maxcvatom,9,"improper:cvatom");
    }
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom && alloc) {
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
  if (cvflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
      cvatom[i][9] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

void Improper::ev_tally(int i1, int i2, int i3, int i4,
                        int nlocal, int newton_bond,
                        double eimproper, double *f1, double *f3, double *f4,
                        double vb1x, double vb1y, double vb1z,
                        double vb2x, double vb2y, double vb2z,
                        double vb3x, double vb3y, double vb3z)
{
  double eimproperquarter,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy += eimproper;
      else {
        eimproperquarter = 0.25*eimproper;
        if (i1 < nlocal) energy += eimproperquarter;
        if (i2 < nlocal) energy += eimproperquarter;
        if (i3 < nlocal) energy += eimproperquarter;
        if (i4 < nlocal) energy += eimproperquarter;
      }
    }
    if (eflag_atom) {
      eimproperquarter = 0.25*eimproper;
      if (newton_bond || i1 < nlocal) eatom[i1] += eimproperquarter;
      if (newton_bond || i2 < nlocal) eatom[i2] += eimproperquarter;
      if (newton_bond || i3 < nlocal) eatom[i3] += eimproperquarter;
      if (newton_bond || i4 < nlocal) eatom[i4] += eimproperquarter;
    }
  }

  if (vflag_either) {
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i1 < nlocal) {
          virial[0] += 0.25*v[0];
          virial[1] += 0.25*v[1];
          virial[2] += 0.25*v[2];
          virial[3] += 0.25*v[3];
          virial[4] += 0.25*v[4];
          virial[5] += 0.25*v[5];
        }
        if (i2 < nlocal) {
          virial[0] += 0.25*v[0];
          virial[1] += 0.25*v[1];
          virial[2] += 0.25*v[2];
          virial[3] += 0.25*v[3];
          virial[4] += 0.25*v[4];
          virial[5] += 0.25*v[5];
        }
        if (i3 < nlocal) {
          virial[0] += 0.25*v[0];
          virial[1] += 0.25*v[1];
          virial[2] += 0.25*v[2];
          virial[3] += 0.25*v[3];
          virial[4] += 0.25*v[4];
          virial[5] += 0.25*v[5];
        }
        if (i4 < nlocal) {
          virial[0] += 0.25*v[0];
          virial[1] += 0.25*v[1];
          virial[2] += 0.25*v[2];
          virial[3] += 0.25*v[3];
          virial[4] += 0.25*v[4];
          virial[5] += 0.25*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        vatom[i1][0] += 0.25*v[0];
        vatom[i1][1] += 0.25*v[1];
        vatom[i1][2] += 0.25*v[2];
        vatom[i1][3] += 0.25*v[3];
        vatom[i1][4] += 0.25*v[4];
        vatom[i1][5] += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
        vatom[i2][0] += 0.25*v[0];
        vatom[i2][1] += 0.25*v[1];
        vatom[i2][2] += 0.25*v[2];
        vatom[i2][3] += 0.25*v[3];
        vatom[i2][4] += 0.25*v[4];
        vatom[i2][5] += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
        vatom[i3][0] += 0.25*v[0];
        vatom[i3][1] += 0.25*v[1];
        vatom[i3][2] += 0.25*v[2];
        vatom[i3][3] += 0.25*v[3];
        vatom[i3][4] += 0.25*v[4];
        vatom[i3][5] += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
        vatom[i4][0] += 0.25*v[0];
        vatom[i4][1] += 0.25*v[1];
        vatom[i4][2] += 0.25*v[2];
        vatom[i4][3] += 0.25*v[3];
        vatom[i4][4] += 0.25*v[4];
        vatom[i4][5] += 0.25*v[5];
      }
    }
  }

  // per-atom centroid virial
  if (cvflag_atom) {

    // r0 = (r1+r2+r3+r4)/4
    // rij = ri-rj
    // total virial = r10*f1 + r20*f2 + r30*f3 + r40*f4
    // vb1: r12
    // vb2: r32
    // vb3: r43

    if (newton_bond || i1 < nlocal) {
      double a1[3];

      // a1 = r10 = (3*r12 - 2*r32 -   r43)/4
      a1[0] = 0.25*(3*vb1x - 2*vb2x - vb3x);
      a1[1] = 0.25*(3*vb1y - 2*vb2y - vb3y);
      a1[2] = 0.25*(3*vb1z - 2*vb2z - vb3z);

      cvatom[i1][0] += a1[0]*f1[0];
      cvatom[i1][1] += a1[1]*f1[1];
      cvatom[i1][2] += a1[2]*f1[2];
      cvatom[i1][3] += a1[0]*f1[1];
      cvatom[i1][4] += a1[0]*f1[2];
      cvatom[i1][5] += a1[1]*f1[2];
      cvatom[i1][6] += a1[1]*f1[0];
      cvatom[i1][7] += a1[2]*f1[0];
      cvatom[i1][8] += a1[2]*f1[1];
    }
    if (newton_bond || i2 < nlocal) {
      double a2[3];
      double f2[3];

      // a2 = r20 = ( -r12 - 2*r32 -   r43)/4
      a2[0] = 0.25*(-vb1x - 2*vb2x - vb3x);
      a2[1] = 0.25*(-vb1y - 2*vb2y - vb3y);
      a2[2] = 0.25*(-vb1z - 2*vb2z - vb3z);

      f2[0] = - f1[0] - f3[0] - f4[0];
      f2[1] = - f1[1] - f3[1] - f4[1];
      f2[2] = - f1[2] - f3[2] - f4[2];

      cvatom[i2][0] += a2[0]*f2[0];
      cvatom[i2][1] += a2[1]*f2[1];
      cvatom[i2][2] += a2[2]*f2[2];
      cvatom[i2][3] += a2[0]*f2[1];
      cvatom[i2][4] += a2[0]*f2[2];
      cvatom[i2][5] += a2[1]*f2[2];
      cvatom[i2][6] += a2[1]*f2[0];
      cvatom[i2][7] += a2[2]*f2[0];
      cvatom[i2][8] += a2[2]*f2[1];
    }
    if (newton_bond || i3 < nlocal) {
      double a3[3];

      // a3 = r30 = ( -r12 + 2*r32 -   r43)/4
      a3[0] = 0.25*(-vb1x + 2*vb2x - vb3x);
      a3[1] = 0.25*(-vb1y + 2*vb2y - vb3y);
      a3[2] = 0.25*(-vb1z + 2*vb2z - vb3z);

      cvatom[i3][0] += a3[0]*f3[0];
      cvatom[i3][1] += a3[1]*f3[1];
      cvatom[i3][2] += a3[2]*f3[2];
      cvatom[i3][3] += a3[0]*f3[1];
      cvatom[i3][4] += a3[0]*f3[2];
      cvatom[i3][5] += a3[1]*f3[2];
      cvatom[i3][6] += a3[1]*f3[0];
      cvatom[i3][7] += a3[2]*f3[0];
      cvatom[i3][8] += a3[2]*f3[1];
    }
    if (newton_bond || i4 < nlocal) {
      double a4[3];

      // a4 = r40 = ( -r12 + 2*r32 + 3*r43)/4
      a4[0] = 0.25*(-vb1x + 2*vb2x + 3*vb3x);
      a4[1] = 0.25*(-vb1y + 2*vb2y + 3*vb3y);
      a4[2] = 0.25*(-vb1z + 2*vb2z + 3*vb3z);

      cvatom[i4][0] += a4[0]*f4[0];
      cvatom[i4][1] += a4[1]*f4[1];
      cvatom[i4][2] += a4[2]*f4[2];
      cvatom[i4][3] += a4[0]*f4[1];
      cvatom[i4][4] += a4[0]*f4[2];
      cvatom[i4][5] += a4[1]*f4[2];
      cvatom[i4][6] += a4[1]*f4[0];
      cvatom[i4][7] += a4[2]*f4[0];
      cvatom[i4][8] += a4[2]*f4[1];
    }
  }
}

/* ---------------------------------------------------------------------- */


void Improper::problem(const char *filename, int lineno,
                       int i1, int i2, int i3, int i4)
{
  const auto x = atom->x;
  auto warn = fmt::format("Improper problem: {} {} {} {} {} {}\n",
                          comm->me, update->ntimestep, atom->tag[i1],
                          atom->tag[i2], atom->tag[i3], atom->tag[i4]);
  warn += fmt::format("WARNING:   1st atom: {} {:.8} {:.8} {:.8}\n",
                      comm->me,x[i1][0],x[i1][1],x[i1][2]);
  warn += fmt::format("WARNING:   2nd atom: {} {:.8} {:.8} {:.8}\n",
                      comm->me,x[i2][0],x[i2][1],x[i2][2]);
  warn += fmt::format("WARNING:   3rd atom: {} {:.8} {:.8} {:.8}\n",
                      comm->me,x[i3][0],x[i3][1],x[i3][2]);
  warn += fmt::format("WARNING:   4th atom: {} {:.8} {:.8} {:.8}",
                      comm->me,x[i4][0],x[i4][1],x[i4][2]);
  error->warning(filename, lineno, warn);
}

/* ---------------------------------------------------------------------- */

double Improper::memory_usage()
{
  double bytes = (double)comm->nthreads*maxeatom * sizeof(double);
  bytes += (double)comm->nthreads*maxvatom*6 * sizeof(double);
  return bytes;
}
