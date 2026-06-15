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
#include "constants_oxdna.h"
#include "nucleotide_oxdna.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_oxdna_lrf.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "update.h"

#include "math_extra.h"
#include "pair.h"

#include <cmath>
#include <exception>

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
void BondOxdnaFene::compute_backbone_site(double e1[3], double /*e2*/[3], double /*e3*/[3],
                                          double rbk[3]) const
{
  NucleotideOxdna1 oxdna1;
  oxdna1.backbone_site(e1, nullptr, nullptr, rbk);
}

/* ----------------------------------------------------------------------
   compute function for oxDNA FENE-bond interaction
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */
void BondOxdnaFene::compute(int eflag, int vflag)
{
  int a, b, btemp, in, type;
  int a3ptype, atype, btype, b5ptype;    // tetramer types
  double delf[3], delta[3], deltb[3];    // force, torque increment
  double delr_bkbk[3], ebond, fbond;
  double rsq_bkbk, Deltasq, rlogarg;
  double r_bkbk, rr0, rr0sq;
  // vectors COM-backbone site in lab frame
  double ra_cbk[3], rb_cbk[3];
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

  tagint *id3p = atom->id3p;
  tagint *id5p = atom->id5p;
  int *atomtype = atom->type;

  const double rlogarg_min = 0.2;
  ebond = 0.0;
  ev_init(eflag, vflag);

  // nxyz_xtrct = extracted local unit vectors in lab frame from fix OXDNA/LRF
  nxyz_xtrct = fix_lrf->array_atom;

  // loop over FENE bonds

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][1];
    b = bondlist[in][0];
    type = bondlist[in][2];

    // directionality test: a -> b is 3' -> 5'
    if (atom->tag[b] != id5p[a]) {

      btemp = b;
      b = a;
      a = btemp;
    }

    // a now in 3' direction, b in 5' direction

    ax[0] = nxyz_xtrct[a][0];
    ax[1] = nxyz_xtrct[a][1];
    ax[2] = nxyz_xtrct[a][2];
    ay[0] = nxyz_xtrct[a][3];
    ay[1] = nxyz_xtrct[a][4];
    ay[2] = nxyz_xtrct[a][5];
    az[0] = nxyz_xtrct[a][6];
    az[1] = nxyz_xtrct[a][7];
    az[2] = nxyz_xtrct[a][8];
    bx[0] = nxyz_xtrct[b][0];
    bx[1] = nxyz_xtrct[b][1];
    bx[2] = nxyz_xtrct[b][2];
    by[0] = nxyz_xtrct[b][3];
    by[1] = nxyz_xtrct[b][4];
    by[2] = nxyz_xtrct[b][5];
    bz[0] = nxyz_xtrct[b][6];
    bz[1] = nxyz_xtrct[b][7];
    bz[2] = nxyz_xtrct[b][8];

    // determine tetramer types
    // 3'neighbor a - a - b - 5'neighbor b

    if (id3p[a] != -1) {
      a3ptype = atomtype[atom->map(id3p[a])];
    } else
      a3ptype = 0;

    atype = atomtype[a];
    btype = atomtype[b];

    if (id5p[b] != -1) {
      b5ptype = atomtype[atom->map(id5p[b])];
    } else
      b5ptype = 0;

    // vector COM-backbone site a and b
    compute_backbone_site(ax, ay, az, ra_cbk);
    compute_backbone_site(bx, by, bz, rb_cbk);

    // vector backbone site b to a
    delr_bkbk[0] = x[a][0] + ra_cbk[0] - x[b][0] - rb_cbk[0];
    delr_bkbk[1] = x[a][1] + ra_cbk[1] - x[b][1] - rb_cbk[1];
    delr_bkbk[2] = x[a][2] + ra_cbk[2] - x[b][2] - rb_cbk[2];
    rsq_bkbk = delr_bkbk[0] * delr_bkbk[0] + delr_bkbk[1] * delr_bkbk[1] + delr_bkbk[2] * delr_bkbk[2];
    r_bkbk = sqrt(rsq_bkbk);

    rr0 = r_bkbk - r0[type][a3ptype][atype][btype][b5ptype];
    rr0sq = rr0 * rr0;
    Deltasq = Delta[type][a3ptype][atype][btype][b5ptype] * Delta[type][a3ptype][atype][btype][b5ptype];
    rlogarg = 1.0 - rr0sq / Deltasq;

    // energy
    if (eflag) { ebond = -0.5 * k[type] * log(rlogarg); }

    // switching to capped force for r-r0 -> Delta at
    // r > r_max = r0 + Delta*sqrt(1-rlogarg) OR
    // r < r_min = r0 - Delta*sqrt(1-rlogarg)
    if (rlogarg < rlogarg_min) {
      // issue warning, reset rlogarg and rr0 to cap force
      error->warning(FLERR, "FENE bond too long: {} {} {} {}", update->ntimestep, atom->tag[a], atom->tag[b], r_bkbk);
      rlogarg = rlogarg_min;

      // if overstretched F(r)=F(r_max)=F_max, E(r)=E(r_max)+F_max*(r-r_max)
      if (r_bkbk > r0[type][a3ptype][atype][btype][b5ptype]) {
        rr0 = Delta[type][a3ptype][atype][btype][b5ptype] * sqrt(1.0 - rlogarg);
        // energy
        if (eflag) {
          ebond = -0.5 * k[type] * log(rlogarg) + k[type] * sqrt(1.0 - rlogarg) / rlogarg /
                  Delta[type][a3ptype][atype][btype][b5ptype] * (r_bkbk - r0[type][a3ptype][atype][btype][b5ptype] -
                   Delta[type][a3ptype][atype][btype][b5ptype] * sqrt(1.0 - rlogarg));
        }
      }
      // if overcompressed F(r)=F(r_min)=F_max, E(r)=E(r_min)+F_max*(r_min-r)
      else if (r_bkbk < r0[type][a3ptype][atype][btype][b5ptype]) {
        rr0 = -Delta[type][a3ptype][atype][btype][b5ptype] * sqrt(1.0 - rlogarg);
        // energy
        if (eflag) {
          ebond = -0.5 * k[type] * log(rlogarg) +
              k[type] * sqrt(1.0 - rlogarg) / rlogarg / Delta[type][a3ptype][atype][btype][b5ptype] *
                  (r0[type][a3ptype][atype][btype][b5ptype] -
                   Delta[type][a3ptype][atype][btype][b5ptype] * sqrt(1.0 - rlogarg) - r_bkbk);
        }
      }
    }

    fbond = -k[type] * rr0 / rlogarg / Deltasq / r_bkbk;
    delf[0] = delr_bkbk[0] * fbond;
    delf[1] = delr_bkbk[1] * fbond;
    delf[2] = delr_bkbk[2] * fbond;

    // apply force and torque to each of 2 atoms

    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cbk, delf, delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];
    }

    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cbk, delf, deltb);

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
  int m = atom->ntypes;

  memory->create(k, n + 1, "bond:k");
  memory->create(Delta, n + 1, m + 1, m + 1, m + 1, m + 1, "bond:Delta");
  memory->create(r0, n + 1, m + 1, m + 1, m + 1, m + 1, "bond:r0");
  memory->create(setflag, n + 1, "bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondOxdnaFene::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene, oxdna2/fene or oxrna2/fene" + utils::errorurl(21));

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one;
  double Delta_one;
  double r0_one;

  if (narg == 4) {
    k_one = utils::numeric(FLERR, arg[1], false, lmp);
    Delta_one = utils::numeric(FLERR, arg[2], false, lmp);
    r0_one = utils::numeric(FLERR, arg[3], false, lmp);
  } else {
    if (comm->me == 0) {    // read values from potential file
      PotentialFileReader reader(lmp, arg[1], "oxdna potential", " (fene)");
      char *line;
      std::string iloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && potential_name == "fene") {
            k_one = values.next_double();
            Delta_one = values.next_double();
            r0_one = values.next_double();

            break;
          } else
            continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxdna, oxdna2 or oxrna2 potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (potential_name != "fene"))
        error->one(FLERR, "No corresponding fene potential found in file {} for bond type {}", arg[1], arg[0]);
    }

    MPI_Bcast(&k_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&Delta_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&r0_one, 1, MPI_DOUBLE, 0, world);
  }

  int count = 0;
  int n = atom->ntypes;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    for (int n1 = 0; n1 <= n; n1++) {    // type 0 for terminal n2
      for (int n2 = 0; n2 <= n; n2++) {
        for (int n3 = 0; n3 <= n; n3++) {
          for (int n4 = 0; n4 <= n; n4++) {    // type 0 for terminal n3
            Delta[i][n1][n2][n3][n4] = Delta_one;
            r0[i][n1][n2][n3][n4] = r0_one;
          }
        }
      }
    }
    setflag[i] = 1;
    count++;
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene, oxdna2/fene or oxrna2/fene" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   set special_bond settings and check if valid
------------------------------------------------------------------------- */

void BondOxdnaFene::init_style()
{
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR, "Must use 'special_bonds lj 0 1 1' with bond style oxdna/fene, oxdna2/fene, " "oxdna3/fene or oxrna2/fene");

  fix_lrf = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF");
  if (fixes.size() == 0)
    error->all(FLERR, "Fix OXDNA/LRF not found. Ensure pair oxdna/excv is present");
  else
    fix_lrf = dynamic_cast<FixOxdnaLRF *>(fixes[0]);
}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::equilibrium_distance(int i)
{
  return r0[i][0][0][0][0];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondOxdnaFene::write_restart(FILE *fp)
{
  int ii, jj, kk, ll;

  fwrite(&k[0], sizeof(double), 1, fp);
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  for (ii = 1; ii <= atom->ntypes; ii++) {
    for (jj = 1; jj <= atom->ntypes; jj++) {
      for (kk = 1; kk <= atom->ntypes; kk++) {
        for (ll = 1; ll <= atom->ntypes; ll++) {
          fwrite(&Delta[0][ii][jj][kk][ll], sizeof(double), 1, fp);
          fwrite(&Delta[1][ii][jj][kk][ll], sizeof(double), atom->nbondtypes, fp);
          fwrite(&r0[0][ii][jj][kk][ll], sizeof(double), 1, fp);
          fwrite(&r0[1][ii][jj][kk][ll], sizeof(double), atom->nbondtypes, fp);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondOxdnaFene::read_restart(FILE *fp)
{
  int ii, jj, kk, ll;

  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[0], sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    for (ii = 1; ii <= atom->ntypes; ii++) {
      for (jj = 1; jj <= atom->ntypes; jj++) {
        for (kk = 1; kk <= atom->ntypes; kk++) {
          for (ll = 1; ll <= atom->ntypes; ll++) {
            utils::sfread(FLERR, &Delta[0][ii][jj][kk][ll], sizeof(double), 1, fp, nullptr, error);
            utils::sfread(FLERR, &Delta[1][ii][jj][kk][ll], sizeof(double), atom->nbondtypes, fp,
                          nullptr, error);
            utils::sfread(FLERR, &r0[0][ii][jj][kk][ll], sizeof(double), 1, fp, nullptr, error);
            utils::sfread(FLERR, &r0[1][ii][jj][kk][ll], sizeof(double), atom->nbondtypes, fp,
                          nullptr, error);
          }
        }
      }
    }
  }

  MPI_Bcast(&k[0], 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  for (ii = 1; ii <= atom->ntypes; ii++) {
    for (jj = 1; jj <= atom->ntypes; jj++) {
      for (kk = 1; kk <= atom->ntypes; kk++) {
        for (ll = 1; ll <= atom->ntypes; ll++) {
          MPI_Bcast(&Delta[0][ii][jj][kk][ll], 1, MPI_DOUBLE, 0, world);
          MPI_Bcast(&Delta[1][ii][jj][kk][ll], atom->nbondtypes, MPI_DOUBLE, 0, world);
          MPI_Bcast(&r0[0][ii][jj][kk][ll], 1, MPI_DOUBLE, 0, world);
          MPI_Bcast(&r0[1][ii][jj][kk][ll], atom->nbondtypes, MPI_DOUBLE, 0, world);
        }
      }
    }
  }

  for (ii = 1; ii <= atom->nbondtypes; ii++) setflag[ii] = 1;
}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r_bkbk = sqrt(rsq);
  double rr0 = r_bkbk - r0[type][0][0][0][0];
  double rr0sq = rr0 * rr0;
  double Deltasq = Delta[type][0][0][0][0] * Delta[type][0][0][0][0];
  double rlogarg = 1.0 - rr0sq / Deltasq;

  // if r -> Delta, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*Delta something serious is wrong, abort

  if (rlogarg < 0.1) {
    error->warning(FLERR, "FENE bond too long: {} {:.8}", update->ntimestep, sqrt(rsq));
    rlogarg = 0.1;
  }

  double eng = -0.5 * k[type] * log(rlogarg);
  fforce = -k[type] * rr0 / rlogarg / Deltasq / r_bkbk;

  return eng;
}
