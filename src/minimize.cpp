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

#include "minimize.h"

#include "citeme.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "min.h"
#include "timer.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Minimize::Minimize(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void Minimize::command(int narg, char **arg)
{
  if (narg != 4)
    error->all(FLERR, "Illegal minimize command: expected 4 arguments but found {}", narg);

  if (domain->box_exist == 0)
    error->all(FLERR, "Minimize command before simulation box is defined");

  // ignore minimize command, if walltime limit was already reached
  if (timer->is_timeout()) return;

  update->etol = utils::numeric(FLERR, arg[0], false, lmp);
  update->ftol = utils::numeric(FLERR, arg[1], false, lmp);
  update->nsteps = utils::inumeric(FLERR, arg[2], false, lmp);
  update->max_eval = utils::inumeric(FLERR, arg[3], false, lmp);

  if (update->etol < 0.0) error->all(FLERR, "Illegal minimize energy tolerance: {}", update->etol);
  if (update->ftol < 0.0) error->all(FLERR, "Illegal minimize force tolerance: {}", update->ftol);

  if (lmp->citeme) lmp->citeme->flush();
  update->whichflag = 2;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  if (update->laststep < 0) error->all(FLERR, "Too many iterations");

  lmp->init();
  timer->init_timeout();
  update->minimize->setup();

  timer->init();
  timer->barrier_start();
  update->minimize->run(update->nsteps);
  timer->barrier_stop();

  update->minimize->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}
