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

/* ------------------------------------------------------------------------
   Contributing author: AUTHOR_NAME_TBD (AFFILIATION_TBD)

   Assign spin masses and create a random distribution of spin velocities at
   a given spin temperature, for use with inertial spin dynamics.  This is
   the spin analogue of the velocity create command; the ordinary velocity
   command is unaffected and still sets the atom velocities.
------------------------------------------------------------------------- */

#include "velocity_tspin.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "random_park.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr double SPIN_EPSILON = 1.0e-4;

/* ---------------------------------------------------------------------- */

VelocityTSpin::VelocityTSpin(LAMMPS *lmp) :
    Command(lmp), igroup(0), groupbit(0), momentum_flag(1), spinmass_flag(0), spinmass(1.0)
{
}

/* ----------------------------------------------------------------------
   velocity/tspin group-ID create T seed keyword value ...
------------------------------------------------------------------------- */

void VelocityTSpin::command(int narg, char **arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "velocity/tspin", error);

  if (domain->box_exist == 0)
    error->all(FLERR, "Velocity/tspin command before simulation box is defined" +
               utils::errorurl(33));
  if (atom->natoms == 0)
    error->all(FLERR, "Velocity/tspin command on a system without atoms");

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR, "Could not find velocity/tspin group ID {}", arg[0]);
  groupbit = group->bitmask[igroup];

  if (strcmp(arg[1], "create") != 0)
    error->all(FLERR, 1, "Unknown velocity/tspin style: {}", arg[1]);

  const double t_desired = utils::numeric(FLERR, arg[2], false, lmp);
  if (t_desired < 0.0) error->all(FLERR, 2, "Velocity/tspin temperature must be >= 0");
  const int seed = utils::inumeric(FLERR, arg[3], false, lmp);
  if (seed <= 0) error->all(FLERR, 3, "Velocity/tspin seed must be > 0");

  momentum_flag = 1;
  spinmass_flag = 0;
  spinmass = 1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "spinmass") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "velocity/tspin spinmass", error);
      spinmass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (spinmass <= 0.0) error->all(FLERR, iarg + 1, "Velocity/tspin spinmass must be > 0");
      spinmass_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "mom") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "velocity/tspin mom", error);
      momentum_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown velocity/tspin keyword: {}", arg[iarg]);
    }
  }

  create(t_desired, seed);
}

/* ----------------------------------------------------------------------
   assign spin masses and draw spin velocities at temperature t_desired
   each atom type is rescaled to t_desired separately
------------------------------------------------------------------------- */

void VelocityTSpin::create(double t_desired, int seed)
{
  if (!atom->tsp_flag) error->all(FLERR, "Velocity/tspin requires atom style tspin");

  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;

  int *local_dof, *global_dof;
  memory->create(local_dof, ntypes + 1, "velocity/tspin:local_dof");
  memory->create(global_dof, ntypes + 1, "velocity/tspin:global_dof");
  for (int itype = 0; itype <= ntypes; itype++) local_dof[itype] = 0;

  // zero all spin velocities, assign spin masses, count spin degrees of freedom

  for (int i = 0; i < nlocal; i++) {
    s_dot[i][0] = s_dot[i][1] = s_dot[i][2] = 0.0;
    if (!(mask[i] & groupbit)) continue;
    if (spinmass_flag) s_mass[i] = spinmass * (rmass ? rmass[i] : mass[type[i]]);
    if (sp[i][3] > SPIN_EPSILON) local_dof[type[i]] += 3;
  }

  MPI_Allreduce(local_dof, global_dof, ntypes + 1, MPI_INT, MPI_SUM, world);

  RanPark random(lmp, seed + comm->me);
  for (int i = 0; i < 100; i++) random.uniform();

  for (int itype = 1; itype <= ntypes; itype++) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit) || (type[i] != itype)) continue;
      if (sp[i][3] <= SPIN_EPSILON) continue;
      s_dot[i][0] = random.gaussian();
      s_dot[i][1] = random.gaussian();
      s_dot[i][2] = random.gaussian();
    }

    if (momentum_flag) zero_mean(itype);

    // rescale to the requested spin temperature

    double local_ke = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit) || (type[i] != itype)) continue;
      if (sp[i][3] <= SPIN_EPSILON) continue;
      local_ke += s_mass[i] *
          (s_dot[i][0] * s_dot[i][0] + s_dot[i][1] * s_dot[i][1] + s_dot[i][2] * s_dot[i][2]);
    }

    double global_ke = 0.0;
    MPI_Allreduce(&local_ke, &global_ke, 1, MPI_DOUBLE, MPI_SUM, world);

    if (global_dof[itype] == 0) continue;
    const double cur_t = global_ke * force->mvv2e / (global_dof[itype] * force->boltz);
    if (cur_t <= 0.0) continue;

    const double rescale = sqrt(t_desired / cur_t);
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit) || (type[i] != itype)) continue;
      if (sp[i][3] <= SPIN_EPSILON) continue;
      s_dot[i][0] *= rescale;
      s_dot[i][1] *= rescale;
      s_dot[i][2] *= rescale;
    }
  }

  memory->destroy(local_dof);
  memory->destroy(global_dof);
}

/* ----------------------------------------------------------------------
   subtract the mean spin velocity of atom type itype
------------------------------------------------------------------------- */

void VelocityTSpin::zero_mean(int itype)
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  int *type = atom->type;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double local_sum[4] = {0.0, 0.0, 0.0, 0.0};

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || (type[i] != itype)) continue;
    if (sp[i][3] <= SPIN_EPSILON) continue;
    local_sum[0] += s_dot[i][0];
    local_sum[1] += s_dot[i][1];
    local_sum[2] += s_dot[i][2];
    local_sum[3] += 1.0;
  }

  double global_sum[4];
  MPI_Allreduce(local_sum, global_sum, 4, MPI_DOUBLE, MPI_SUM, world);
  if (global_sum[3] == 0.0) return;

  const double invn = 1.0 / global_sum[3];
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || (type[i] != itype)) continue;
    if (sp[i][3] <= SPIN_EPSILON) continue;
    s_dot[i][0] -= global_sum[0] * invn;
    s_dot[i][1] -= global_sum[1] * invn;
    s_dot[i][2] -= global_sum[2] * invn;
  }
}
