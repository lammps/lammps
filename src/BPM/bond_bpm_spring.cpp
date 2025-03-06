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

#include "bond_bpm_spring.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBPMSpring::BondBPMSpring(LAMMPS *_lmp) :
    BondBPM(_lmp), k(nullptr), av(nullptr), ecrit(nullptr), gamma(nullptr),
    id_fix_property_bond(nullptr), vol_current(nullptr), dvol0(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  volume_flag = 0;
  writedata = 0;

  nhistory = 1;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_SPRING");

  single_extra = 1;
  svector = new double[1];

  nmax = 0;

  comm_forward = 0;
  comm_reverse = 0;
}

/* ---------------------------------------------------------------------- */

BondBPMSpring::~BondBPMSpring()
{
  delete[] svector;
  if (id_fix_property_bond && modify->nfix) {
    modify->delete_fix(id_fix_property_bond);
    delete[] id_fix_property_bond;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(ecrit);
    memory->destroy(gamma);
    memory->destroy(av);
  }

  memory->destroy(dvol0);
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMSpring::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r;
  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  r = sqrt(delx * delx + dely * dely + delz * delz);
  bondstore[n][0] = r;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) { fix_bond_history->update_atom_value(i, m, 0, r); }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) { fix_bond_history->update_atom_value(j, m, 0, r); }
    }
  }

  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMSpring::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      if (type < 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);
    }
  }

  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::compute(int eflag, int vflag)
{
  int i, bond_change_flag;
  double *vol0, *vol;

  if (volume_flag) {
    vol0 = atom->dvector[index_vol0];
    vol = atom->dvector[index_vol];

    // grow + initialize dvol0 as necessary
    if (nmax < atom->nmax) {
      nmax = atom->nmax;
      memory->create(dvol0, nmax, "bond/bpm/spring:dvol0");
      for (i = 0; i < nmax; i++) dvol0[i] = 0.0;
    }
  }

  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();

    if (volume_flag) {
      vol_current = vol0;
      bond_change_flag = calculate_vol();

      // zero dvol0, not needed since vol0 just calculated
      for (i = 0; i < nmax; i++) dvol0[i] = 0.0;
    }
  }

  if (volume_flag) {
    vol_current = vol;
    bond_change_flag = calculate_vol();

    // Update vol0 to account for any new bonds
    if (bond_change_flag) {
      update_vol0();

      // forward new vol0 to ghosts before force calculation
      vol_current = vol0;
      comm->forward_comm(this);

      bond_change_flag = 0;
    }
  }

  if (hybrid_flag) fix_bond_history->compress_history();

  int i1, i2, itmp, n, type;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, rsq, r, r0, rinv, smooth, fbond, dot;
  double vol_sum, vol0_sum, vol_temp;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dim = domain->dimension;
  double invdim = 1.0 / dim;

  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero
    if (r0 < EPSILON || std::isnan(r0)) r0 = store_bond(n, i1, i2);

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    e = (r - r0) / r0;

    if ((fabs(e) > ecrit[type]) && break_flag) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);

      if (volume_flag) {
        bond_change_flag = 1;
        vol_temp = r0 * r0;
        if (dim == 3) vol_temp *= r0;
        if (newton_bond || i1 < nlocal) dvol0[i1] -= vol_temp;
        if (newton_bond || i2 < nlocal) dvol0[i2] -= vol_temp;
      }

      continue;
    }

    rinv = 1.0 / r;
    if (normalize_flag)
      fbond = -k[type] * e;
    else
      fbond = k[type] * (r0 - r);

    if (volume_flag) {
      vol_sum = vol[i1] + vol[i2];
      vol0_sum = vol0[i1] + vol0[i2];
      fbond += av[type] * (pow(vol_sum / vol0_sum, invdim) - 1.0 - e);
    }

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

  // Update vol0 to account for any broken bonds
  if (volume_flag && bond_change_flag) update_vol0();

  if (hybrid_flag) fix_bond_history->uncompress_history();
}

/* ---------------------------------------------------------------------- */

int BondBPMSpring::calculate_vol()
{
  int n, i1, i2;
  double r0, delx, dely, delz, rsq, vol_temp;

  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;
  int newton_bond = force->newton_bond;
  int dim = domain->dimension;

  double **x = atom->x;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < ntotal; n++) vol_current[n] = 0.0;

  int bond_change_flag = 0;

  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    r0 = bondstore[n][0];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    rsq = delx * delx + dely * dely + delz * delz;

    vol_temp = rsq;
    if (dim == 3) vol_temp *= sqrt(rsq);

    if (newton_bond || i1 < nlocal) vol_current[i1] += vol_temp;
    if (newton_bond || i2 < nlocal) vol_current[i2] += vol_temp;

    // If bond hasn't been set - increment dvol0 too to update vol0
    if (r0 < EPSILON || std::isnan(r0)) {
      bond_change_flag = 1;
      if (newton_bond || i1 < nlocal) dvol0[i1] += vol_temp;
      if (newton_bond || i2 < nlocal) dvol0[i2] += vol_temp;
    }
  }

  if (newton_bond) comm->reverse_comm(this);
  comm->forward_comm(this);

  return bond_change_flag;
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::update_vol0()
{
  // accumulate changes in vol0 from ghosts
  vol_current = dvol0;
  if (force->newton_bond) comm->reverse_comm(this);

  double *vol0 = atom->dvector[index_vol0];
  for (int i = 0; i < atom->nlocal; i++) vol0[i] += dvol0[i];

  // zero dvol0 for next change
  for (int i = 0; i < nmax; i++) dvol0[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");
  memory->create(ecrit, np1, "bond:ecrit");
  memory->create(gamma, np1, "bond:gamma");
  memory->create(av, np1, "bond:av");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMSpring::coeff(int narg, char **arg)
{
  if ((!volume_flag && narg != 4) || (volume_flag && narg != 5))
    error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double ecrit_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);

  double av_one = 0.0;
  if (volume_flag) av_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    av[i] = av_one;
    setflag[i] = 1;
    count++;

    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMSpring::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/spring requires ghost atoms store velocity");

  if (volume_flag && !id_fix_property_bond) {
    id_fix_property_bond = utils::strdup("BOND_BPM_SPRING_FIX_PROPERTY_ATOM");
    modify->add_fix(fmt::format("{} all property/atom d_vol d_vol0 ghost yes writedata no",
                                id_fix_property_bond));

    int tmp1 = 0, tmp2 = 0;
    index_vol = atom->find_custom("vol", tmp1, tmp2);
    index_vol0 = atom->find_custom("vol0", tmp1, tmp2);
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::settings(int narg, char **arg)
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
    } else if (strcmp(arg[iarg], "volume/factor") == 0) {
      if (iarg + 1 > narg)
        error->all(FLERR, "Illegal bond bpm command, missing option for volume/factor");
      volume_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);

      if (volume_flag) {
        comm_forward = 1;
        comm_reverse = 1;
      } else {
        comm_forward = 0;
        comm_reverse = 0;
      }

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

void BondBPMSpring::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gamma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&av[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMSpring::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &av[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&av[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMSpring::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
  fwrite(&volume_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMSpring::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &volume_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&volume_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMSpring::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  double r0;
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) r0 = fix_bond_history->get_atom_value(i, n, 0);
  }

  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double e = (r - r0) / r0;

  if (normalize_flag)
    fforce = -k[type] * e;
  else
    fforce = k[type] * (r0 - r);

  if (volume_flag) {
    double invdim = 1.0 / domain->dimension;
    double *vol0 = atom->dvector[index_vol0];
    double *vol = atom->dvector[index_vol];
    double vol_sum = vol[i] + vol[j];
    double vol0_sum = vol0[i] + vol0[j];
    fforce += av[type] * (pow(vol_sum / vol0_sum, invdim) - 1.0 - e);
  }

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
    double smooth = (r - r0) / (r0 * ecrit[type]);
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = r0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int BondBPMSpring::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = vol_current[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    vol_current[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int BondBPMSpring::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = vol_current[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) vol_current[i] = buf[m++];
}
