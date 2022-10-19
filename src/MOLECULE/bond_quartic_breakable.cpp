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

/* ----------------------------------------------------------------------
   Contributing authors: Chris Lorenz and Mark Stevens (SNL)
                         Tiedong Sun (NTU)
------------------------------------------------------------------------- */

#include "bond_quartic_breakable.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_CUBEROOT2;

/* ---------------------------------------------------------------------- */

BondQuarticBreakable::BondQuarticBreakable(LAMMPS *_lmp) :
    Bond(_lmp), k(nullptr), b1(nullptr), b2(nullptr), rc(nullptr), u0(nullptr)
{
  partial_flag = 1;
  breakable_flag = 1;
}

/* ---------------------------------------------------------------------- */

BondQuarticBreakable::~BondQuarticBreakable()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(b1);
    memory->destroy(b2);
    memory->destroy(rc);
    memory->destroy(u0);
  }
}

/* ---------------------------------------------------------------------- */

void BondQuarticBreakable::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (vflag_global == VIRIAL_FDOTR) force->pair->vflag_either = force->pair->vflag_global = 1;

  if (breakable_flag) {
    if (evflag) {
      if (eflag) eval<1,1,1>();
      else eval<1,0,1>();
    } else eval<1,0,0>();
  } else {
    if (evflag) {
      if (eflag) eval<0,1,1>();
      else eval<0,0,1>();
    } else eval<0,0,0>();
  }
}

/* ---------------------------------------------------------------------- */
template <int BREAKABLE, int EFLAG, int EVFLAG>
void BondQuarticBreakable::eval()
{
  int i1, i2, n, m, type, itype, jtype;
  double delx, dely, delz, ebond, fbond, evdwl, fpair;
  double r, rsq, dr, r2, ra, rb, sr2, sr6;

  ebond = evdwl = sr6 = 0.0;

  double **cutsq = force->pair->cutsq;
  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken

    if (BREAKABLE && bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;

    // quartic term
    if (rsq > rc[type] * rc[type]) {
      // for breakable bonds
      // if bond breaks, set type to 0
      //   both in temporary bondlist and permanent bond_type
      // if this proc owns both atoms,
      //   negate bond_type twice if other atom stores it
      // if other proc owns 2nd atom, other proc will also break bond
      if (BREAKABLE) {
        bondlist[n][2] = 0;
        for (m = 0; m < atom->num_bond[i1]; m++)
          if (atom->bond_atom[i1][m] == atom->tag[i2]) atom->bond_type[i1][m] = 0;
        if (i2 < atom->nlocal)
          for (m = 0; m < atom->num_bond[i2]; m++)
            if (atom->bond_atom[i2][m] == atom->tag[i1]) atom->bond_type[i2][m] = 0;
        fbond = 0.0;
        if (EFLAG) ebond = 0.0;
      } else {
        // unbreakable bond
        fbond = 0.0;
        if (EFLAG) ebond = u0[type];
      }
    } else {
        r = sqrt(rsq);
        dr = r - rc[type];
        r2 = dr * dr;
        ra = dr - b1[type];
        rb = dr - b2[type];
        fbond = -k[type] / r * (r2 * (ra + rb) + 2.0 * dr * ra * rb);
        if (EFLAG) ebond = k[type] * r2 * ra * rb + u0[type];
    }

    // LJ term
    if (rsq < MY_CUBEROOT2) {
      sr2 = 1.0 / rsq;
      sr6 = sr2 * sr2 * sr2;
      fbond += 48.0 * sr6 * (sr6 - 0.5) / rsq;
      if (EFLAG) ebond += 4.0 * sr6 * (sr6 - 1.0) + 1.0;
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    if (EVFLAG) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);

    // subtract out pairwise contribution from 2 atoms via pair->single()
    // required since special_bond = 1,1,1
    // tally energy/virial in pair, using newton_bond as newton flag

    if (BREAKABLE) {
      itype = atom->type[i1];
      jtype = atom->type[i2];

      if (rsq < cutsq[itype][jtype]) {
        evdwl = -force->pair->single(i1, i2, itype, jtype, rsq, 1.0, 1.0, fpair);
        fpair = -fpair;

        if (newton_bond || i1 < nlocal) {
          f[i1][0] += delx * fpair;
          f[i1][1] += dely * fpair;
          f[i1][2] += delz * fpair;
        }
        if (newton_bond || i2 < nlocal) {
          f[i2][0] -= delx * fpair;
          f[i2][1] -= dely * fpair;
          f[i2][2] -= delz * fpair;
        }

        if (EVFLAG)
          force->pair->ev_tally(i1, i2, nlocal, newton_bond, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
   }
}

/* ---------------------------------------------------------------------- */

void BondQuarticBreakable::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");
  memory->create(b1, np1, "bond:b1");
  memory->create(b2, np1, "bond:b2");
  memory->create(rc, np1, "bond:rc");
  memory->create(u0, np1, "bond:u0");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondQuarticBreakable::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double b1_one = utils::numeric(FLERR, arg[2], false, lmp);
  double b2_one = utils::numeric(FLERR, arg[3], false, lmp);
  double rc_one = utils::numeric(FLERR, arg[4], false, lmp);
  double u0_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    b1[i] = b1_one;
    b2[i] = b2_one;
    rc[i] = rc_one;
    u0[i] = u0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if pair defined and special_bond settings are valid
------------------------------------------------------------------------- */

void BondQuarticBreakable::init_style()
{
  if (breakable_flag) {
    if (force->pair == nullptr || force->pair->single_enable == 0)
      error->all(FLERR, "Pair style does not support bond_style quartic");
    if (force->angle || force->dihedral || force->improper)
      error->all(FLERR, "Bond style quartic cannot be used with 3,4-body interactions");
    if (atom->molecular == Atom::TEMPLATE)
      error->all(FLERR, "Bond style quartic cannot be used with atom style template");

    // special bonds must be 1 1 1

    if (force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
      error->all(FLERR, "Bond style quartic requires special_bonds = 1,1,1");
  }
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondQuarticBreakable::equilibrium_distance(int i)
{
  // equilibrium distance of the quartic potential, excluding LJ
  double quartic_equ = MY_CUBEROOT2;
  double b2_4ac = 9.0 * (b1[i] * b1[i] + b2[i] *  b2[i]) - 14.0 * b1[i] * b2[i];

  if (b2_4ac >= 0.0) {
    quartic_equ = (8.0 * rc[i] + 3.0 * (b1[i] + b2[i]) - sqrt(b2_4ac) ) / 8.0;
  } else return rc[i] < MY_CUBEROOT2 ? MY_CUBEROOT2 : rc[i];

  // when quartic_equ < MY_CUBEROOT2, the bond length hits the LJ wall, using MY_CUBEROOT2
  // as aproximate bond equilibrium length
  return quartic_equ < MY_CUBEROOT2 ? MY_CUBEROOT2 : quartic_equ;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondQuarticBreakable::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&b1[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&b2[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&rc[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&u0[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondQuarticBreakable::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &b1[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &b2[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &rc[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &u0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&b1[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&b2[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rc[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&u0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondQuarticBreakable::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g %g %g %g %g\n", i, k[i], b1[i], b2[i], rc[i], u0[i]);
}

/* ---------------------------------------------------------------------- */

double BondQuarticBreakable::single(int type, double rsq, int i, int j, double &fforce)
{
  double r, dr, r2, ra, rb, sr2, sr6;

  // broken bond
  if (type <= 0) {
    fforce = 0.0;
    return 0.0;
  }
  // unbreakable bond, r > rc
  if (rsq > rc[type] * rc[type]) {
    fforce = 0.0;
    return u0[type];
  }

  double eng = 0.0;

  // for breakable bond
  // subtract out pairwise contribution from 2 atoms via pair->single()
  // required since special_bond = 1,1,1

  if (breakable_flag) {
    int itype = atom->type[i];
    int jtype = atom->type[j];

    if (rsq < force->pair->cutsq[itype][jtype]) {
      double tmp;
      eng = -force->pair->single(i, j, itype, jtype, rsq, 1.0, 1.0, tmp);
    }
  }

  // quartic bond
  // 1st portion is from quartic term
  // 2nd portion is from LJ term cut at 2^(1/6) with eps = sigma = 1.0

  r = sqrt(rsq);
  dr = r - rc[type];
  r2 = dr * dr;
  ra = dr - b1[type];
  rb = dr - b2[type];

  eng += k[type] * r2 * ra * rb + u0[type];
  fforce = -k[type] / r * (r2 * (ra + rb) + 2.0 * dr * ra * rb);

  if (rsq < MY_CUBEROOT2) {
    sr2 = 1.0 / rsq;
    sr6 = sr2 * sr2 * sr2;
    eng += 4.0 * sr6 * (sr6 - 1.0) + 1.0;
    fforce += 48.0 * sr6 * (sr6 - 0.5) / rsq;
  }

  return eng;
}
