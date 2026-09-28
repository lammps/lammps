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

#include "fix_propel_ou.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropelOU::FixPropelOU(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), decay(0.0), kick(0.0), ilevel_respa(0), factive(nullptr), rng(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix propel/ou", error);
  if (narg > 6) error->all(FLERR, 6, "Unknown fix propel/ou keyword {}", arg[6]);

  magnitude = utils::numeric(FLERR, arg[3], false, lmp);
  if (magnitude <= 0.0) error->all(FLERR, 3, "Fix propel/ou magnitude must be > 0");

  tau = utils::numeric(FLERR, arg[4], false, lmp);
  if (tau <= 0.0) error->all(FLERR, 4, "Fix propel/ou correlation time must be > 0");

  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  if (seed <= 0) error->all(FLERR, 5, "Fix propel/ou seed must be > 0");

  restart_global = 1;
  restart_peratom = 1;
  create_attribute = 1;
  respa_level_support = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_virial = 1;
  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;
  maxexchange = 3;

  // initialize Marsaglia RNG with processor-unique seed

  rng = new RanMars(lmp, seed + comm->me);

  // per-atom active force, registered with the Atom class

  FixPropelOU::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // draw the initial active force of each atom from the stationary distribution

  const int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) FixPropelOU::set_arrays(i);
}

/* ---------------------------------------------------------------------- */

FixPropelOU::~FixPropelOU()
{
  if (copymode) return;

  atom->delete_callback(id, Atom::GROW);
  atom->delete_callback(id, Atom::RESTART);
  memory->destroy(factive);
  delete rng;
}

/* ---------------------------------------------------------------------- */

int FixPropelOU::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropelOU::init()
{
  reset_dt();

  if (utils::strmatch(update->integrate_style, "^respa")) {
    int max_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    ilevel_respa = max_respa;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, max_respa);
  }
}

/* ----------------------------------------------------------------------
   prefactors of the exact update of the Ornstein-Uhlenbeck process over
   one timestep: f(t+dt) = decay * f(t) + kick * gaussian(0,1)
------------------------------------------------------------------------- */

void FixPropelOU::reset_dt()
{
  decay = exp(-update->dt / tau);
  kick = magnitude * sqrt(1.0 - decay * decay);
}

/* ----------------------------------------------------------------------
   apply the current active force without advancing the process, so that
   consecutive runs continue the same trajectory as a single longer run
------------------------------------------------------------------------- */

void FixPropelOU::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    apply_force(vflag);
  } else {
    auto *respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    apply_force(vflag);
    respa->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixPropelOU::post_force(int vflag)
{
  advance();
  apply_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPropelOU::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ----------------------------------------------------------------------
   advance the active force of all atoms in the group by one timestep
------------------------------------------------------------------------- */

void FixPropelOU::advance()
{
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int dim = domain->dimension;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    factive[i][0] = decay * factive[i][0] + kick * rng->gaussian();
    factive[i][1] = decay * factive[i][1] + kick * rng->gaussian();
    if (dim == 3) factive[i][2] = decay * factive[i][2] + kick * rng->gaussian();
  }
}

/* ----------------------------------------------------------------------
   add the active force to the atoms in the group and tally its virial
   contribution using unwrapped coordinates
------------------------------------------------------------------------- */

void FixPropelOU::apply_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  const int nlocal = atom->nlocal;

  double vi[6], unwrap[3];
  if (vflag)
    v_setup(vflag);
  else
    evflag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double fx = factive[i][0];
    const double fy = factive[i][1];
    const double fz = factive[i][2];
    f[i][0] += fx;
    f[i][1] += fy;
    f[i][2] += fz;

    if (evflag) {
      domain->unmap(x[i], image[i], unwrap);
      vi[0] = fx * unwrap[0];
      vi[1] = fy * unwrap[1];
      vi[2] = fz * unwrap[2];
      vi[3] = fx * unwrap[1];
      vi[4] = fx * unwrap[2];
      vi[5] = fy * unwrap[2];
      v_tally(i, vi);
    }
  }
}

/* ----------------------------------------------------------------------
   pack the per-processor RNG state into the restart file so that a run
   continued from a restart reproduces the original stochastic trajectory
------------------------------------------------------------------------- */

void FixPropelOU::write_restart(FILE *fp)
{
  int nsize = RanMars::STATE_SIZE * comm->nprocs + 1;    // pRNG state per proc + nprocs
  auto *list = new double[nsize];

  if (comm->me == 0) list[0] = comm->nprocs;

  double state[RanMars::STATE_SIZE];
  rng->get_state(state);
  MPI_Gather(state, RanMars::STATE_SIZE, MPI_DOUBLE, list + 1, RanMars::STATE_SIZE, MPI_DOUBLE, 0,
             world);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), nsize, fp);
  }
  delete[] list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restore the RNG state
------------------------------------------------------------------------- */

void FixPropelOU::restart(char *buf)
{
  auto *list = (double *) buf;

  int nprocs = (int) list[0];
  if (nprocs != comm->nprocs) {
    if (comm->me == 0)
      error->warning(FLERR, "Different number of procs. Cannot restore RNG state.");
  } else {
    // the size of the stored states depends on the version that wrote the restart file
    const int stride = RanMars::state_size(list + 1);
    rng->set_state(list + 1 + comm->me * stride);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropelOU::memory_usage()
{
  return (double) atom->nmax * 3 * sizeof(double);
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPropelOU::grow_arrays(int nmax)
{
  memory->grow(factive, nmax, 3, "propel/ou:factive");
  array_atom = factive;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixPropelOU::copy_arrays(int i, int j, int /*delflag*/)
{
  factive[j][0] = factive[i][0];
  factive[j][1] = factive[i][1];
  factive[j][2] = factive[i][2];
}

/* ----------------------------------------------------------------------
   initialize the active force of a new atom: draw it from the stationary
   distribution of the process for atoms in the group, else zero it
------------------------------------------------------------------------- */

void FixPropelOU::set_arrays(int i)
{
  if (atom->mask[i] & groupbit) {
    factive[i][0] = magnitude * rng->gaussian();
    factive[i][1] = magnitude * rng->gaussian();
    factive[i][2] = (domain->dimension == 3) ? magnitude * rng->gaussian() : 0.0;
  } else {
    factive[i][0] = factive[i][1] = factive[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropelOU::pack_exchange(int i, double *buf)
{
  buf[0] = factive[i][0];
  buf[1] = factive[i][1];
  buf[2] = factive[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropelOU::unpack_exchange(int nlocal, double *buf)
{
  factive[nlocal][0] = buf[0];
  factive[nlocal][1] = buf[1];
  factive[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropelOU::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it
  buf[0] = 4;
  buf[1] = factive[i][0];
  buf[2] = factive[i][1];
  buf[3] = factive[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropelOU::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  factive[nlocal][0] = extra[nlocal][m++];
  factive[nlocal][1] = extra[nlocal][m++];
  factive[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropelOU::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropelOU::size_restart(int /*nlocal*/)
{
  return 4;
}
