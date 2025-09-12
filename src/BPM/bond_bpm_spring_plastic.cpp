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
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "bond_bpm_spring_plastic.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

static const char cite_bpm_plastic[] =
    "BPM/spring/plastic bond style: doi:10.1016/j.powtec.2024.120563\n\n"
    "@Article{Clemmer2025,\n"
    " author =  {Clemmer, Joel T. and Lechman, Jeremy B.},\n"
    " title =   {Onset and impact of plastic deformation in granular compaction},\n"
    " journal = {Powder Technology},\n"
    " year =    2025,\n"
    " volume =  452,\n"
    " number =  28,\n"
    " pages =   {120563}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

BondBPMSpringPlastic::BondBPMSpringPlastic(LAMMPS *_lmp) :
    BondBPM(_lmp), k(nullptr), eplastic(nullptr), ecrit(nullptr), gamma(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  writedata = 0;

  nhistory = 2;
  update_flag = 1;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_SPRING_PLASTIC");

  single_extra = 2;
  svector = new double[2];

  if (lmp->citeme) lmp->citeme->add(cite_bpm_plastic);
}

/* ---------------------------------------------------------------------- */

BondBPMSpringPlastic::~BondBPMSpringPlastic()
{
  delete[] svector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(ecrit);
    memory->destroy(gamma);
    memory->destroy(eplastic);
  }
}

/* ----------------------------------------------------------------------
  Store data for all bonds, called once
------------------------------------------------------------------------- */

void BondBPMSpringPlastic::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      if (type <= 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(FLERR, delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, 0);
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMSpringPlastic::compute(int eflag, int vflag)
{
  pre_compute();

  int i1, i2, itmp, n, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, ep, rsq, r, r0, rinv, smooth, fbond, dot;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;
  const bool allow_breaks = (update->setupflag == 0) && break_flag;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    // If bond hasn't been set (should be initialized to zero)
    r0 = bondstore[n][0];
    if (r0 < EPSILON || std::isnan(r0)) {
      r0 = bondstore[n][0] = r;
      bondstore[n][1] = 0.0;
      process_new(n, i1, i2);
    }
    ep = bondstore[n][1];

    e = (r - r0) / r0;

    if ((fabs(e) > ecrit[type]) && allow_breaks) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    if (e > (ep + eplastic[type])) {
      ep = e - eplastic[type];
      bondstore[n][1] = ep;
    }

    if (e < (ep - eplastic[type])) {
      ep = e + eplastic[type];
      bondstore[n][1] = ep;
    }

    rinv = 1.0 / r;
    if (normalize_flag)
      fbond = -k[type] * (e - ep);
    else
      fbond = -k[type] * (r - (1.0 + ep) * r0);

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx * delvx + dely * delvy + delz * delvz;
    fbond -= gamma[type] * dot * rinv;
    fbond *= rinv;

    if (smooth_flag) {
      smooth = (r - r0) / (r0 * ecrit[type]);
      smooth *= smooth;
      smooth *= smooth;
      smooth *= smooth;
      smooth = 1 - smooth;
      fbond *= smooth;
    }

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

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, 0.0, fbond, delx, dely, delz);
  }

  post_compute();
}

/* ---------------------------------------------------------------------- */

void BondBPMSpringPlastic::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");
  memory->create(ecrit, np1, "bond:ecrit");
  memory->create(gamma, np1, "bond:gamma");
  memory->create(eplastic, np1, "bond:eplastic");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMSpringPlastic::coeff(int narg, char **arg)
{
  if (narg != 5)
    error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double ecrit_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double eplastic_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    eplastic[i] = eplastic_one;
    setflag[i] = 1;
    count++;

    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMSpringPlastic::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/spring requires ghost atoms store velocity");
}

/* ---------------------------------------------------------------------- */

void BondBPMSpringPlastic::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  int iarg;
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "smooth") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for smooth");
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 > narg)
        error->all(FLERR, "Illegal bond bpm command, missing option for normalize");
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond bpm command, invalid argument {}", arg[iarg]);
    }
  }

  if (smooth_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off smoothing with break no option");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMSpringPlastic::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&eplastic[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMSpringPlastic::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &eplastic[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&eplastic[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMSpringPlastic::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMSpringPlastic::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMSpringPlastic::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  // ep can be updated, so search bondlist vs. fix_bond_history->get_atom_value()
  //   slower than other bpm/bond styles' single method
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **bondstore = fix_bond_history->bondstore;

  tagint tagi = tag[i];
  tagint tagj = tag[j];
  tagint tag1, tag2;
  double r0 = 0.0, ep = 0.0;
  for (int n = 0; n < nbondlist; n++) {
    tag1 = tag[bondlist[n][0]];
    tag2 = tag[bondlist[n][1]];
    if ((tag1 == tagi && tag2 == tagj) || (tag1 == tagj && tag2 == tagi)) {
      r0 = bondstore[n][0];
      ep = bondstore[n][1];
      break;
    }
  }

  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double e = (r0 != 0.0) ? (r - r0) / r0 : 0.0;

  if (normalize_flag)
    fforce = -k[type] * (e - ep);
  else
    fforce = -k[type] * (r - (1.0 + ep) * r0);

  double **x = atom->x;
  double **v = atom->v;
  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];
  double delvx = v[i][0] - v[j][0];
  double delvy = v[i][1] - v[j][1];
  double delvz = v[i][2] - v[j][2];
  double dot = delx * delvx + dely * delvy + delz * delvz;
  fforce -= gamma[type] * dot * rinv;
  fforce *= rinv;

  if (smooth_flag) {
    double smooth = (r0 != 0.0) ? (r - r0) / (r0 * ecrit[type]) : 0.0;
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = r0;
  svector[1] = (1.0 + ep) * r0;

  return 0.0;
}
