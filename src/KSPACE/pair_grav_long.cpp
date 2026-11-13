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
----------------------------------
   Contributing author: Joel Clemmer (SNL)

   Real-space part of the Ewald sum for Newtonian gravity,
   E = -G m_i m_j erfc(g_ewald r) / r, to be used with kspace_style ewald/grav.
   Adapted from pair style coul/long.

   Reference: L. Hernquist, F. R. Bouchet, Y. Suto, ApJS 75, 231 (1991).

   Note: Eq. (2.14b) as printed has 2 alpha / sqrt(pi) in front of the
   exponential, without the factor r = |x - nL|, which must be a misprint
------------------------------------------------------------------------- */

#include "pair_grav_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairGravLong::PairGravLong(LAMMPS *lmp) :
    Pair(lmp), gconst(0.0), cut_grav(0.0), cut_gravsq(0.0), g_ewald(0.0)
{
  ewaldflag = pppmflag = 1;
  ncoultablebits = 0;    // erfc() is evaluated directly, no tables
}

/* ---------------------------------------------------------------------- */

PairGravLong::~PairGravLong()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairGravLong::compute(int eflag, int vflag)
{
  double egrav = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  const double *const rmass = atom->rmass;
  const double *const mass = atom->mass;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *const special_coul = force->special_coul;
  const int newton_pair = force->newton_pair;

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double mtmp = rmass ? rmass[i] : mass[type[i]];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_gravsq) {
        const double mj = rmass ? rmass[j] : mass[type[j]];
        const double r = sqrt(rsq);
        const double grij = g_ewald * r;
        const double expm2 = exp(-grij * grij);
        const double erfcr = erfc(grij);
        const double prefactor = -gconst * mtmp * mj / r;
        double forcegrav = prefactor * (erfcr + EWALD_F * grij * expm2);
        if (factor < 1.0) forcegrav -= (1.0 - factor) * prefactor;

        const double fpair = forcegrav / rsq;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          egrav = prefactor * erfcr;
          if (factor < 1.0) egrav -= (1.0 - factor) * prefactor;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, egrav, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGravLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   gravitational constant G in the current units: 1.0 for lj units, else
   G_SI * mass^2 / (length * energy) with the units of mass, length, and
   energy expressed in SI, from CODATA 2022 values
   https://physics.nist.gov/cuu/Constants/Table/allascii.txt
------------------------------------------------------------------------- */

double PairGravLong::default_gconst() const
{
  constexpr double G_SI = 6.67430e-11;              // m^3 / (kg s^2)
  constexpr double AMU = 1.66053906892e-27;         // kg, also g/mol per particle
  constexpr double NA = 6.02214076e23;              // 1/mol
  constexpr double KCAL_MOL = 4184.0 / NA;          // J
  constexpr double EV = 1.602176634e-19;            // J
  constexpr double HARTREE = 4.3597447222060e-18;   // J
  constexpr double BOHR = 5.29177210544e-11;        // m

  const char *style = update->unit_style;
  if (strcmp(style, "lj") == 0) return 1.0;
  if (strcmp(style, "real") == 0) return G_SI * AMU * AMU / (1.0e-10 * KCAL_MOL);
  if (strcmp(style, "metal") == 0) return G_SI * AMU * AMU / (1.0e-10 * EV);
  if (strcmp(style, "si") == 0) return G_SI;
  if (strcmp(style, "cgs") == 0) return G_SI * 1.0e-3 * 1.0e-3 / (1.0e-2 * 1.0e-7);
  if (strcmp(style, "electron") == 0) return G_SI * AMU * AMU / (BOHR * HARTREE);
  if (strcmp(style, "micro") == 0) return G_SI * 1.0e-15 * 1.0e-15 / (1.0e-6 * 1.0e-15);
  if (strcmp(style, "nano") == 0) return G_SI * 1.0e-21 * 1.0e-21 / (1.0e-9 * 1.0e-21);
  error->all(FLERR, "Pair style grav/long has no gravitational constant for units {}; "
             "set it with the gconst keyword", style);
  return 0.0;
}

/* ----------------------------------------------------------------------
   global settings: cutoff and optional gravitational constant
------------------------------------------------------------------------- */

void PairGravLong::settings(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "pair_style grav/long", error);

  cut_grav = utils::numeric(FLERR, arg[0], false, lmp);
  if (cut_grav <= 0.0) error->all(FLERR, 1, "Cutoff for pair_style grav/long must be > 0");

  int gflag = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "gconst") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style grav/long gconst", error);
      gconst = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (gconst <= 0.0)
        error->all(FLERR, iarg + 1,
                   "Gravitational constant for pair_style grav/long must be > 0");
      gflag = 1;
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown pair_style grav/long keyword: {}", arg[iarg]);
    }
  }
  if (!gflag) gconst = default_gconst();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGravLong::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGravLong::init_style()
{
  neighbor->add_request(this);

  cut_gravsq = cut_grav * cut_grav;

  // ensure use of KSpace long-range solver, set g_ewald

  if ((force->kspace == nullptr) || !utils::strmatch(force->kspace_style, "/grav"))
    error->all(FLERR, "Pair style grav/long requires a gravity KSpace style like ewald/grav");
  g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGravLong::init_one(int /*i*/, int /*j*/)
{
  return cut_grav;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGravLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) fwrite(&setflag[i][j], sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGravLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGravLong::write_restart_settings(FILE *fp)
{
  fwrite(&gconst, sizeof(double), 1, fp);
  fwrite(&cut_grav, sizeof(double), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGravLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &gconst, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_grav, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&gconst, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_grav, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairGravLong::single(int i, int j, int itype, int jtype, double rsq, double factor_coul,
                            double /*factor_lj*/, double &fforce)
{
  const double mi = atom->rmass ? atom->rmass[i] : atom->mass[itype];
  const double mj = atom->rmass ? atom->rmass[j] : atom->mass[jtype];
  const double r = sqrt(rsq);
  const double grij = g_ewald * r;
  const double expm2 = exp(-grij * grij);
  const double erfcr = erfc(grij);
  const double prefactor = -gconst * mi * mj / r;

  double forcegrav = prefactor * (erfcr + EWALD_F * grij * expm2);
  double phigrav = prefactor * erfcr;
  if (factor_coul < 1.0) {
    forcegrav -= (1.0 - factor_coul) * prefactor;
    phigrav -= (1.0 - factor_coul) * prefactor;
  }
  fforce = forcegrav / rsq;
  return phigrav;
}

/* ----------------------------------------------------------------------
   "cut_coul" is the name under which the KSpace styles query the
   real-space cutoff of the Ewald splitting
------------------------------------------------------------------------- */

void *PairGravLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "cut_coul") == 0) return (void *) &cut_grav;
  if (strcmp(str, "G_grav") == 0) return (void *) &gconst;
  return nullptr;
}
