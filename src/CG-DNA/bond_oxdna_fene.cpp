/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "bond_oxdna_fene.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include "math_extra.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondOxdnaFene::~BondOxdnaFene()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(Delta);
    memory->destroy(r0);
  }
}

/* ----------------------------------------------------------------------
    compute vector COM-sugar-phosphate backbone interaction site in oxDNA
------------------------------------------------------------------------- */
void BondOxdnaFene::compute_interaction_sites(double e1[3], double /*e2*/[3], double /*e3*/[3],
                                              double r[3]) const
{
  constexpr double d_cs = -0.4;

  r[0] = d_cs * e1[0];
  r[1] = d_cs * e1[1];
  r[2] = d_cs * e1[2];
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

void BondOxdnaFene::ev_tally_xyz(int i, int j, int nlocal, int newton_bond, double ebond, double fx,
                                 double fy, double fz, double delx, double dely, double delz)
{
  double ebondhalf, v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond)
        energy += ebond;
      else {
        ebondhalf = 0.5 * ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5 * ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx * fx;
    v[1] = dely * fy;
    v[2] = delz * fz;
    v[3] = delx * fy;
    v[4] = delx * fz;
    v[5] = dely * fz;

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
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5 * v[0];
        vatom[i][1] += 0.5 * v[1];
        vatom[i][2] += 0.5 * v[2];
        vatom[i][3] += 0.5 * v[3];
        vatom[i][4] += 0.5 * v[4];
        vatom[i][5] += 0.5 * v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5 * v[0];
        vatom[j][1] += 0.5 * v[1];
        vatom[j][2] += 0.5 * v[2];
        vatom[j][3] += 0.5 * v[3];
        vatom[j][4] += 0.5 * v[4];
        vatom[j][5] += 0.5 * v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute function for oxDNA FENE-bond interaction
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */
void BondOxdnaFene::compute(int eflag, int vflag)
{
  int a, b, in, type;
  double delf[3], delta[3], deltb[3];    // force, torque increment
  double delr[3], ebond, fbond;
  double rsq, Deltasq, rlogarg;
  double r, rr0, rr0sq;
  // vectors COM-backbone site in lab frame
  double ra_cs[3], rb_cs[3];
  // Cartesian unit vectors in lab frame
  double ax[3], ay[3], az[3];
  double bx[3], by[3], bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  ebond = 0.0;
  ev_init(eflag, vflag);

  // n(x/y/z)_xtrct = extracted local unit vectors in lab frame from oxdna_excv
  int dim;
  nx_xtrct = (double **) force->pair->extract("nx", dim);
  ny_xtrct = (double **) force->pair->extract("ny", dim);
  nz_xtrct = (double **) force->pair->extract("nz", dim);

  // loop over FENE bonds

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][1];
    b = bondlist[in][0];
    type = bondlist[in][2];

    ax[0] = nx_xtrct[a][0];
    ax[1] = nx_xtrct[a][1];
    ax[2] = nx_xtrct[a][2];
    ay[0] = ny_xtrct[a][0];
    ay[1] = ny_xtrct[a][1];
    ay[2] = ny_xtrct[a][2];
    az[0] = nz_xtrct[a][0];
    az[1] = nz_xtrct[a][1];
    az[2] = nz_xtrct[a][2];
    bx[0] = nx_xtrct[b][0];
    bx[1] = nx_xtrct[b][1];
    bx[2] = nx_xtrct[b][2];
    by[0] = ny_xtrct[b][0];
    by[1] = ny_xtrct[b][1];
    by[2] = ny_xtrct[b][2];
    bz[0] = nz_xtrct[b][0];
    bz[1] = nz_xtrct[b][1];
    bz[2] = nz_xtrct[b][2];

    // vector COM-backbone site a and b
    compute_interaction_sites(ax, ay, az, ra_cs);
    compute_interaction_sites(bx, by, bz, rb_cs);

    // vector backbone site b to a
    delr[0] = x[a][0] + ra_cs[0] - x[b][0] - rb_cs[0];
    delr[1] = x[a][1] + ra_cs[1] - x[b][1] - rb_cs[1];
    delr[2] = x[a][2] + ra_cs[2] - x[b][2] - rb_cs[2];
    rsq = delr[0] * delr[0] + delr[1] * delr[1] + delr[2] * delr[2];
    r = sqrt(rsq);

    rr0 = r - r0[type];
    rr0sq = rr0 * rr0;
    Deltasq = Delta[type] * Delta[type];
    rlogarg = 1.0 - rr0sq / Deltasq;

    // if r -> Delta, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*Delta something serious is wrong, abort

    if (rlogarg < 0.1) {
      error->warning(FLERR, "FENE bond too long: {} {} {} {}", update->ntimestep, atom->tag[a],
                     atom->tag[b], r);
      rlogarg = 0.1;
    }

    fbond = -k[type] * rr0 / rlogarg / Deltasq / r;
    delf[0] = delr[0] * fbond;
    delf[1] = delr[1] * fbond;
    delf[2] = delr[2] * fbond;

    // energy

    if (eflag) { ebond = -0.5 * k[type] * log(rlogarg); }

    // apply force and torque to each of 2 atoms

    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cs, delf, delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];
    }

    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cs, delf, deltb);

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];
    }

    // increment energy and virial
    // NOTE: The virial is calculated on the 'molecular' basis.
    // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

    if (evflag)
      ev_tally_xyz(a, b, nlocal, newton_bond, ebond, delf[0], delf[1], delf[2], x[a][0] - x[b][0],
                   x[a][1] - x[b][1], x[a][2] - x[b][2]);
  }
}

/* ---------------------------------------------------------------------- */

void BondOxdnaFene::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k, n + 1, "bond:k");
  memory->create(Delta, n + 1, "bond:Delta");
  memory->create(r0, n + 1, "bond:r0");
  memory->create(setflag, n + 1, "bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondOxdnaFene::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double Delta_one = utils::numeric(FLERR, arg[2], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[3], false, lmp);

  int count = 0;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    Delta[i] = Delta_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene");
}

/* ----------------------------------------------------------------------
   set special_bond settings and check if valid
------------------------------------------------------------------------- */

void BondOxdnaFene::init_style()
{
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(
        FLERR,
        "Must use 'special_bonds lj 0 1 1' with bond style oxdna/fene, oxdna2/fene or oxrna2/fene");
}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondOxdnaFene::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Delta[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondOxdnaFene::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Delta[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Delta[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondOxdnaFene::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, k[i], r0[i], Delta[i]);
}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r = sqrt(rsq);
  double rr0 = r - r0[type];
  double rr0sq = rr0 * rr0;
  double Deltasq = Delta[type] * Delta[type];
  double rlogarg = 1.0 - rr0sq / Deltasq;

  // if r -> Delta, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*Delta something serious is wrong, abort

  if (rlogarg < 0.1) {
    error->warning(FLERR, "FENE bond too long: {} {:.8}", update->ntimestep, sqrt(rsq));
    rlogarg = 0.1;
  }

  double eng = -0.5 * k[type] * log(rlogarg);
  fforce = -k[type] * rr0 / rlogarg / Deltasq / r;

  return eng;
}
