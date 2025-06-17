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

#include "fix_set.h"

#include "atom.h"
#include "error.h"
#include "set.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SETCOMMAND,FIXSET};     // also used in Set class

/* ---------------------------------------------------------------------- */

FixSet::FixSet(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR, 1, "Illegal fix set command: need at least eight arguments");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Fix {} Nevery must be > 0", style);

  reneighbor = utils::inumeric(FLERR, arg[4], false, lmp);
  if (reneighbor < 0 || reneighbor > 1)
    error->all(FLERR, "Fix {} rnflag must be 0/1", style);

  // create instance of Set class

  set = new Set(lmp);

  // pass remaining args to Set class
  // only keywords which use per-atom variables are currently allowed
  // NOTE: could also allow when set style = region,
  //       since atoms may move in/out of regions

  set->process_args(FIXSET,narg-5,&arg[5]);
}

/* ---------------------------------------------------------------------- */

FixSet::~FixSet()
{
  delete set;
}

/* ---------------------------------------------------------------------- */

int FixSet::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   use the Set instance to update per-atom properties
   NOTE: could return count of updated atoms from Set for use as fix output
---------------------------------------------------------------------- */

void FixSet::end_of_step()
{
  // select which atoms to act on

  set->selection(atom->nlocal);

  // loop over list of actions to reset atom attributes

  set->invoke_actions();

  // trigger reneighboring on next timestep if requested

  if (reneighbor) next_reneighbor = update->ntimestep + 1;
}
