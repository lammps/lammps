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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "bond_rheo_shell.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_surface.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "fix_rheo.h"
#include "fix_rheo_oxidation.h"
#include "fix_store_local.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;
using namespace RHEO_NS;

/* ---------------------------------------------------------------------- */

BondRHEOShell::BondRHEOShell(LAMMPS *_lmp) :
    BondBPM(_lmp), k(nullptr), ecrit(nullptr), gamma(nullptr), dbond(nullptr), nbond(nullptr),
    id_fix(nullptr), compute_surface(nullptr)
{
  partial_flag = 1;
  comm_reverse = 1;

  nhistory = 2;
  update_flag = 1;
  id_fix_bond_history = utils::strdup("HISTORY_RHEO_SHELL");
  ignore_special_flag = 1;

  tform = -1;

  single_extra = 1;
  svector = new double[1];

  // For nbond, create an instance of fix property atom
  // Need restarts + exchanging with neighbors since it needs to persist
  // between timesteps (fix property atom will handle callbacks)

  int tmp1, tmp2;
  index_nb = atom->find_custom("shell_nbond", tmp1, tmp2);
  if (index_nb == -1) {
    id_fix = utils::strdup("bond_rheo_shell_fix_property_atom");
    modify->add_fix(fmt::format("{} all property/atom i_shell_nbond", id_fix));
    index_nb = atom->find_custom("shell_nbond", tmp1, tmp2);
  }
  nbond = atom->ivector[index_nb];

  //Store non-persistent per atom quantities, intermediate

  nmax_store = atom->nmax;
  memory->create(dbond, nmax_store, "rheo/react:dbond");
}

/* ---------------------------------------------------------------------- */

BondRHEOShell::~BondRHEOShell()
{
  if (modify->nfix) modify->delete_fix(id_fix);
  delete[] id_fix;
  delete[] svector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(ecrit);
    memory->destroy(gamma);
  }

  memory->destroy(dbond);
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondRHEOShell::store_bond(int n, int i, int j)
{
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  bondstore[n][0] = 0.0;
  bondstore[n][1] = 0.0;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) {
        fix_bond_history->update_atom_value(i, m, 0, 0.0);
        fix_bond_history->update_atom_value(i, m, 1, 0.0);
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) {
        fix_bond_history->update_atom_value(j, m, 0, 0.0);
        fix_bond_history->update_atom_value(j, m, 1, 0.0);
      }
    }
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondRHEOShell::store_data()
{
  int i, j, m, type;
  int **bond_type = atom->bond_type;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      if (type < 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      fix_bond_history->update_atom_value(i, m, 0, 0.0);
      fix_bond_history->update_atom_value(i, m, 1, 0.0);
    }
  }

  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondRHEOShell::compute(int eflag, int vflag)
{
  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, rsq, r, r0, rinv, dr, fbond, dot, t;
  double dt = update->dt;

  ev_init(eflag, vflag);

  double *rsurface = compute_surface->rsurface;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *status = atom->rheo_status;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;

  if (atom->nmax > nmax_store) {
    nmax_store = atom->nmax;
    memory->destroy(dbond);
    memory->create(dbond, nmax_store, "rheo/shell:dbond");
  }

  size_t nbytes = nmax_store * sizeof(int);
  memset(&dbond[0], 0, nbytes);

  for (n = 0; n < nbondlist; n++) {
    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];
    t = bondstore[n][1];

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

    // If bond hasn't been set - zero data
    if (t < EPSILON || std::isnan(t)) t = store_bond(n, i1, i2);

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    // Bond has not yet formed, check if in range + update timer
    if (t < tform) {

      // Check if eligible
      if (r > rmax || rsurface[i1] > rsurf || rsurface[i2] > rsurf) {
        bondlist[n][2] = 0;
        process_ineligibility(i1, i2);
        continue;
      }

      // Check ellapsed time
      t += dt;
      bondstore[n][1] = t;
      if (t >= tform) {
        bondstore[n][0] = r;
        r0 = r;
        if (newton_bond || i1 < nlocal) dbond[i1]++;
        if (newton_bond || i2 < nlocal) dbond[i2]++;
      } else {
        continue;
      }
    }

    e = (r - r0) / r0;
    if (fabs(e) > ecrit[type]) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      if (newton_bond || i1 < nlocal) dbond[i1]--;
      if (newton_bond || i2 < nlocal) dbond[i2]--;
      continue;
    }

    rinv = 1.0 / r;
    dr = r - r0;
    fbond = 2 * k[type] * (-dr + dr * dr * dr / (r0 * r0 * ecrit[type] * ecrit[type]));

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx * delvx + dely * delvy + delz * delvz;
    fbond -= gamma[type] * dot * rinv;
    fbond *= rinv;

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

  // Communicate changes in nbond
  if (newton_bond) comm->reverse_comm(this);

  for (int i = 0; i < nlocal; i++) {
    nbond[i] += dbond[i];

    // If it has bonds, no shifting
    if (nbond[i] != 0) status[i] |= STATUS_NO_SHIFT;
  }

  if (hybrid_flag) fix_bond_history->uncompress_history();
}

/* ---------------------------------------------------------------------- */

void BondRHEOShell::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");
  memory->create(ecrit, np1, "bond:ecrit");
  memory->create(gamma, np1, "bond:gamma");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondRHEOShell::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double ecrit_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    setflag[i] = 1;
    count++;

    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondRHEOShell::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond rheo/shell requires ghost atoms store velocity");

  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use bond rheo/shell");
  class FixRHEO *fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (!fix_rheo->surface_flag)
    error->all(FLERR, "Bond rheo/shell requires surface calculation in fix rheo");
  compute_surface = fix_rheo->compute_surface;

  fixes = modify->get_fix_by_style("^rheo/oxidation$");
  if (fixes.size() == 0)
    error->all(FLERR, "Need to define fix rheo/oxidation to use bond rheo/shell");
  class FixRHEOOxidation *fix_rheo_oxidation = dynamic_cast<FixRHEOOxidation *>(fixes[0]);

  rsurf = fix_rheo_oxidation->rsurf;
  rmax = fix_rheo_oxidation->cut;
}

/* ---------------------------------------------------------------------- */

void BondRHEOShell::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  int iarg;
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "t/form") == 0) {
      if (iarg + 1 > narg) utils::missing_cmd_args(FLERR, "bond rheo/shell t/form", error);
      tform = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond rheo/shell command, invalid argument {}", arg[iarg]);
    }
  }

  if (tform < 0.0)
    error->all(FLERR, "Illegal bond rheo/shell command, must specify positive formation time");
}

/* ----------------------------------------------------------------------
   used to check bond communiction cutoff - not perfect, estimates based on local-local only
------------------------------------------------------------------------- */

double BondRHEOShell::equilibrium_distance(int /*i*/)
{
  // Divide out heuristic prefactor added in comm class
  return max_stretch * rmax / 1.5;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondRHEOShell::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondRHEOShell::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondRHEOShell::write_restart_settings(FILE *fp)
{
  fwrite(&tform, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondRHEOShell::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) { utils::sfread(FLERR, &tform, sizeof(double), 1, fp, nullptr, error); }
  MPI_Bcast(&tform, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

int BondRHEOShell::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;

  for (i = first; i < last; i++) { buf[m++] = dbond[i]; }
  return m;
}

/* ---------------------------------------------------------------------- */

void BondRHEOShell::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dbond[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double BondRHEOShell::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  double r0, t;
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) {
      r0 = fix_bond_history->get_atom_value(i, n, 0);
      t = fix_bond_history->get_atom_value(i, n, 1);
    }
  }

  svector[1] = t;
  if (t < tform) return 0.0;

  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double dr = r0 - r;
  fforce = 2 * k[type] * (dr + dr * dr * dr / (r0 * r0 * ecrit[type] * ecrit[type]));

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

  // set single_extra quantities

  svector[0] = r0;

  return 0.0;
}

/* ----------------------------------------------------------------------
    Similar to BondBPM->process_broken(), but don't send to FixStoreLocal
 ------------------------------------------------------------------------- */

void BondRHEOShell::process_ineligibility(int i, int j)
{
  // Manually search and remove from atom arrays
  int m, n;
  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  if (i < nlocal) {
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == tag[j] && setflag[bond_type[i][m]]) {
        bond_type[i][m] = 0;
        n = num_bond[i];
        bond_type[i][m] = bond_type[i][n - 1];
        bond_atom[i][m] = bond_atom[i][n - 1];
        for (auto &ihistory : histories) {
          auto fix_bond_history2 = dynamic_cast<FixBondHistory *>(ihistory);
          fix_bond_history2->shift_history(i, m, n - 1);
          fix_bond_history2->delete_history(i, n - 1);
        }
        num_bond[i]--;
        break;
      }
    }
  }

  if (j < nlocal) {
    for (m = 0; m < num_bond[j]; m++) {
      if (bond_atom[j][m] == tag[i] && setflag[bond_type[j][m]]) {
        bond_type[j][m] = 0;
        n = num_bond[j];
        bond_type[j][m] = bond_type[j][n - 1];
        bond_atom[j][m] = bond_atom[j][n - 1];
        for (auto &ihistory : histories) {
          auto fix_bond_history2 = dynamic_cast<FixBondHistory *>(ihistory);
          fix_bond_history2->shift_history(j, m, n - 1);
          fix_bond_history2->delete_history(j, n - 1);
        }
        num_bond[j]--;
        break;
      }
    }
  }
}
