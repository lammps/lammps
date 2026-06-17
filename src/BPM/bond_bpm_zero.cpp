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

#include "bond_bpm_zero.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBPMZero::BondBPMZero(LAMMPS *_lmp) :
    BondBPM(_lmp), ecrit(nullptr), id_fix_property_bond(nullptr)
{
  partial_flag = 1;
  writedata = 0;

  // Whether a manybody term is needed, depends on specific model
  //   Here toggled by a placeholder setting to demonstrate an arbitrary calculation
  manybody_flag = 0;

  update_flag = 0; // Whether stored values can evolve
  nhistory = 1; // Number of values stored per bond
  id_fix_bond_history = utils::strdup("HISTORY_BPM_ZERO");

  single_extra = 1;
  svector = new double[1];

  nmax = 0;

  comm_forward = 0;
  comm_reverse = 0;
}

/* ---------------------------------------------------------------------- */

BondBPMZero::~BondBPMZero()
{
  delete[] svector;
  if (id_fix_property_bond && modify->nfix) {
    modify->delete_fix(id_fix_property_bond);
    delete[] id_fix_property_bond;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(ecrit);
  }
}

/* ----------------------------------------------------------------------
  Store data for all bonds, called once
------------------------------------------------------------------------- */

void BondBPMZero::store_data()
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
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMZero::compute(int eflag, int vflag)
{
  // Rearrange stored bond values for hybrid bond styles, handled by parent class
  pre_compute();

  double *manybody_term;

  if (manybody_flag) {
    manybody_term = atom->dvector[index_manybody];
    calculate_manybody();
  }

  int i1, i2, itmp, n, type;
  double delx, dely, delz, rsq, r, r0, e, fbond, ebond;

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
      // redefine stored values to current value
      r0 = bondstore[n][0] = r;

      // add new bond and save stored values
      process_new(n, i1, i2);
    }

    // Calculate strain to see if bond breaks
    //   (may use other breakage criterion)
    e = (r - r0) / r0;

    if ((fabs(e) > ecrit[type]) && allow_breaks) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    // Calculate forces and energies
    fbond = 0.0;
    ebond = 0.0;

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

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }

  // Revert changes for hybrid bond style, handled by parent class
  post_compute();
}

/* ----------------------------------------------------------------------
  If needed, calculate manybody term
     here, calculate sum of r0 as an example
------------------------------------------------------------------------- */

void BondBPMZero::calculate_manybody()
{
  int n, i1, i2;
  double r0;

  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;
  int newton_bond = force->newton_bond;
  int dim = domain->dimension;

  double *manybody = atom->dvector[index_manybody];
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < ntotal; n++) manybody[n] = 0.0;

  int bond_change_flag = 0;

  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    r0 = bondstore[n][0];

    if (newton_bond || i1 < nlocal) manybody[i1] += r0;
    if (newton_bond || i2 < nlocal) manybody[i2] += r0;
  }

  if (newton_bond) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void BondBPMZero::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(ecrit, np1, "bond:ecrit");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMZero::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double ecrit_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    ecrit[i] = ecrit_one;
    setflag[i] = 1;
    count++;

    // calculate maximum stretch to check communication cutoffs
    if (1.0 + ecrit[i] > max_stretch) max_stretch = 1.0 + ecrit[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMZero::init_style()
{
  BondBPM::init_style();

  // If using manybody term, create fix property/atom to store result
  if (manybody_flag && !id_fix_property_bond) {
    id_fix_property_bond = utils::strdup("BOND_BPM_ZERO_FIX_PROPERTY_ATOM");
    modify->add_fix(fmt::format("{} all property/atom d_manybody ghost yes writedata no",
                                id_fix_property_bond));

    int tmp1 = 0, tmp2 = 0;
    index_manybody = atom->find_custom("manybody", tmp1, tmp2);
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMZero::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  int iarg;
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "manybody") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "bond_style bpm/zero manybody", error);
      manybody_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);

      // Set communication size to avoid overloading buffers
      if (manybody_flag) {
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
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMZero::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&ecrit[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMZero::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &ecrit[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&ecrit[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMZero::write_restart_settings(FILE *fp)
{
  fwrite(&manybody_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMZero::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &manybody_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&manybody_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMZero::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  double r0;
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) r0 = fix_bond_history->get_atom_value(i, n, 0);
  }

  double r = sqrt(rsq);
  double rinv = 1.0 / r;
  double e = (r - r0) / r0;

  // calculate force and energy

  double ebond = 0.0;
  fforce = 0.0;

  // set single_extra quantities

  svector[0] = r0;

  return ebond;
}

/* ---------------------------------------------------------------------- */

int BondBPMZero::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  double *manybody = atom->dvector[index_manybody];
  int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = manybody[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void BondBPMZero::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  double *manybody = atom->dvector[index_manybody];
  for (int i = 0; i < n; i++) {
    int j = list[i];
    manybody[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int BondBPMZero::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  double *manybody = atom->dvector[index_manybody];
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = manybody[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void BondBPMZero::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  double *manybody = atom->dvector[index_manybody];
  int last = first + n;
  for (int i = first; i < last; i++) {
    manybody[i] = buf[m++];
  }
}
