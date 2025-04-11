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
   Package      FixPIMDBLangevin

   Purpose      Path Integral Molecular Dynamics of Bosons with Langevin Thermostat
   Copyright    Hirshberg lab @ Tel Aviv University
   Authors      Ofir Blumer, Jacob Higer, Yotam Feldman (yotam.feldman at gmail.com), Barak Hirshberg (hirshb at tauex.tau.ac.il)

   Updated      Jan-06-2025
   Version      1.0
------------------------------------------------------------------------- */

#include "fix_pimd_langevin_bosonic.h"

#include "bosonic_exchange.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::FixPIMDBLangevin(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDLangevin(lmp, narg, filtered_args = filter_args(narg, arg)), filtered_narg(narg),
    nbosons(atom->nlocal)
{
  bosonic_exchange = new BosonicExchange(lmp, atom->nlocal, np, universe->me, true, false);
  synch_energies = true;

  // Loop over the arguments with i++ instead of i+=2 like the parent,
  // because the parent's loop has i++ inside some if blocks.
  for (int i = 3; i < narg - 1; i++) {
    if ((strcmp(arg[i], "method") == 0) && (strcmp(arg[i + 1], "pimd") != 0)) {
      error->universe_all(FLERR, "Method not supported in fix pimdb/langevin; only method PIMD");
    } else if (strcmp(arg[i], "scale") == 0) {
      error->universe_all(FLERR,
                          "The scale parameter of the PILE_L thermostat is not supported for "
                          "pimdb, and should be removed.");
    } else if ((strcmp(arg[i], "iso") == 0) || (strcmp(arg[i], "aniso") == 0) ||
               (strcmp(arg[i], "barostat") == 0) || (strcmp(arg[i], "taup") == 0)) {
      error->universe_all(FLERR, "Barostat parameters are not available for pimdb.");
    } else if (strcmp(arg[i], "esynch") == 0) {
      if (strcmp(arg[i + 1], "yes") == 0) {
        synch_energies = true;
      } else if (strcmp(arg[i + 1], "no") == 0) {
        synch_energies = false;
      } else {
        error->universe_all(FLERR, "The esynch parameter can only receive yes or no!");
      }
    }
  }

  if (fmmode != PHYSICAL) {
    error->universe_all(
        FLERR,
        "The only available fmmode for pimdb is physical, please remove the fmmode keyword.");
  }
  if (ensemble != NVE && ensemble != NVT) {
    error->universe_all(FLERR,
                        "The only available ensembles for pimdb are nve and nvt, please choose one "
                        "of these ensembles.");
  }

  method = PIMD;
  size_vector = 6;
  memory->create(f_tag_order, nbosons, 3, "FixPIMDBLangevin:f_tag_order");

  if (cmode != SINGLE_PROC)
    error->universe_all(FLERR,
                        fmt::format("Fix {} only supports a single processor per bead", style));
}

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::~FixPIMDBLangevin()
{
  memory->destroy(f_tag_order);
  for (int i = 0; i < filtered_narg; ++i) delete[] filtered_args[i];
  delete[] filtered_args;
  delete bosonic_exchange;
}

/* ---------------------------------------------------------------------- */

char **FixPIMDBLangevin::filter_args(int narg, char **arg)
{
  filtered_narg = narg;
  char **filtered_args = new char *[narg];
  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i], "esynch") == 0) {
      filtered_args[i] = utils::strdup("");
    } else {
      filtered_args[i] = utils::strdup(arg[i]);
    }
  }
  return filtered_args;
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::prepare_coordinates()
{
  inter_replica_comm(atom->x);
  double ff = fbond * atom->mass[atom->type[0]];
  int nlocal = atom->nlocal;
  double *me_bead_positions = *(atom->x);
  double *last_bead_positions = &bufsortedall[x_last * nlocal][0];
  double *next_bead_positions = &bufsortedall[x_next * nlocal][0];

  bosonic_exchange->prepare_with_coordinates(me_bead_positions, last_bead_positions,
                                             next_bead_positions, beta_np, ff);
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::spring_force()
{

  for (int i = 0; i < nbosons; i++) {
    f_tag_order[i][0] = 0.0;
    f_tag_order[i][1] = 0.0;
    f_tag_order[i][2] = 0.0;
  }
  bosonic_exchange->spring_force(f_tag_order);

  double **f = atom->f;
  tagint *tag = atom->tag;
  for (int i = 0; i < nbosons; i++) {
    f[i][0] += f_tag_order[tag[i] - 1][0];
    f[i][1] += f_tag_order[tag[i] - 1][1];
    f[i][2] += f_tag_order[tag[i] - 1][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::compute_spring_energy()
{
  se_bead = bosonic_exchange->get_bead_spring_energy();

  if (synch_energies) {
    MPI_Allreduce(&se_bead, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    total_spring_energy /= universe->procs_per_world[universe->iworld];
  } else {
    total_spring_energy = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::compute_t_prim()
{
  if (synch_energies) {
    double prim = bosonic_exchange->prim_estimator();
    MPI_Allreduce(&prim, &t_prim, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  } else {
    t_prim = bosonic_exchange->prim_estimator();
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMDBLangevin::compute_vector(int n)
{
  if (0 <= n && n < 6) {
    return FixPIMDLangevin::compute_vector(n);
  } else {
    error->universe_all(FLERR, "Fix only has 6 outputs!");
  }
  return 0.0;
}
