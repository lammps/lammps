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

#include "fix_pimd_nvt_bosonic.h"

#include "bosonic_exchange.h"

#include "atom.h"
#include "error.h"
#include "memory.h"
#include "universe.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;

namespace {
enum { PHYSICAL, NORMAL };
}

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::FixPIMDBNVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg, true), bosonic_exchange(nullptr), f_tag_order(nullptr),
    nbosons(atom->nlocal)
{
  method = PIMD;
  mtchain = 2;
  t_period = 100.0;

  parse_nvt_arguments(narg, arg, [this](int parse_narg, char **parse_arg, int &i) {
    return parse_bosonic_keyword(parse_narg, parse_arg, i);
  });

  if (method != PIMD)
    error->all(FLERR, "Fix pimd/nvt/bosonic only supports method pimd");
  if (fmmode != PHYSICAL)
    error->all(FLERR, "Fix pimd/nvt/bosonic only supports fmmode physical");

  finish_nuclear_constructor_setup();

  bosonic_exchange = new BosonicExchange(lmp, nbosons, np, universe->me, true, false);
  memory->create(f_tag_order, nbosons, 3, "FixPIMDBNVT:f_tag_order");
}

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::~FixPIMDBNVT()
{
  memory->destroy(f_tag_order);
  delete bosonic_exchange;
}

/* ---------------------------------------------------------------------- */

bool FixPIMDBNVT::parse_bosonic_keyword(int narg, char **arg, int &i)
{
  if (strcmp(arg[i], "nhc") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, "fix pimd/nvt/bosonic nhc", error);
    mtchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
    if (mtchain < 1) error->all(FLERR, "Nose-Hoover chain length for fix {} must be >= 1", style);
    i += 2;
    return true;
  }
  return false;
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::spring_force()
{
  double *me_bead_positions = *(atom->x);
  double *last_bead_positions = &bufsortedall[x_last * nbosons][0];
  double *next_bead_positions = &bufsortedall[x_next * nbosons][0];
  double ff = fbond * atom->mass[atom->type[0]];

  bosonic_exchange->prepare_with_coordinates(me_bead_positions, last_bead_positions,
                                             next_bead_positions, beta_np, ff);

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

void FixPIMDBNVT::compute_spring_energy()
{
  se_bead = bosonic_exchange->get_bead_spring_energy();
  MPI_Allreduce(&se_bead, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::compute_t_prim()
{
  double prim = bosonic_exchange->prim_estimator();
  MPI_Allreduce(&prim, &t_prim, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}
