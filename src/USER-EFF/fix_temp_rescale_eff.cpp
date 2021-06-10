// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero
------------------------------------------------------------------------- */

#include "fix_temp_rescale_eff.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixTempRescaleEff::FixTempRescaleEff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix temp/rescale/eff command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix temp/rescale/eff command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;
  ecouple_flag = 1;

  t_start = utils::numeric(FLERR,arg[4],false,lmp);
  t_stop = utils::numeric(FLERR,arg[5],false,lmp);
  t_window = utils::numeric(FLERR,arg[6],false,lmp);
  fraction = utils::numeric(FLERR,arg[7],false,lmp);

  // create a new compute temp/eff
  // id = fix-ID + temp, compute group = fix group

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp/eff",
                                  id_temp,group->names[igroup]));
  tflag = 1;

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempRescaleEff::~FixTempRescaleEff()
{
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempRescaleEff::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleEff::init()
{
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/rescale/eff does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleEff::end_of_step()
{
  double t_current = temperature->compute_scalar();
  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/rescale/eff cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int *spin = atom->spin;
    double *ervel = atom->ervel;

    if (which == NOBIAS) {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          v[i][0] *= factor;
          v[i][1] *= factor;
          v[i][2] *= factor;
          if (abs(spin[i])==1)
            ervel[i] *= factor;
        }
      }
    } else {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          temperature->remove_bias(i,v[i]);
          v[i][0] *= factor;
          v[i][1] *= factor;
          v[i][2] *= factor;
          if (abs(spin[i])==1)
            ervel[i] *= factor;
          temperature->restore_bias(i,v[i]);
        }
      }
    }

  }
}

/* ---------------------------------------------------------------------- */

int FixTempRescaleEff::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp/eff") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleEff::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempRescaleEff::compute_scalar()
{
  return energy;
}
