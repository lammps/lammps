/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_rescale.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "temperature.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

FixTempRescale::FixTempRescale(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 8) error->all("Illegal fix temp/rescale command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix temp/rescale command");

  t_start = atof(arg[4]);
  t_end = atof(arg[5]);
  t_window = atof(arg[6]);
  fraction = atof(arg[7]);

  // optional args

  iregion = -1;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      for (iregion = 0; iregion < domain->nregion; iregion++)
	if (strcmp(arg[iarg+1],domain->regions[iregion]->id) == 0) break;
      if (iregion == domain->nregion)
	error->all("Fix temp/rescale region ID does not exist");
      iarg += 2;
    } else error->all("Illegal fix temp/rescale command");
  }

  // create a new temperature full or region style with fix ID and fix group

  char **newarg = new char*[4];
  newarg[0] = id;
  newarg[1] = group->names[igroup];
  if (iregion == -1) {
    newarg[2] = "full";
    force->add_temp(3,newarg,1);
  } else {
    newarg[2] = "region";
    newarg[3] = domain->regions[iregion]->id;
    force->add_temp(4,newarg,1);
  }
  delete [] newarg;

  temperature = force->find_temp(id);
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::init()
{
  energy = 0.0;
  efactor = (0.5 * force->boltz * temperature->dof);
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::end_of_step()
{
  double t_current = temperature->compute();
  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_end-t_start);

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);

    double **x = atom->x;
    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    if (iregion == -1) {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
	}
      }

    } else {
      efactor = (0.5 * force->boltz * temperature->dof);
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit &&
	    domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
	}
      }
    }
  } else energy = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    temperature = force->find_temp(arg[1]);
    if (temperature == NULL)
      error->all("Could not find fix_modify temperature ID");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning("Group for fix_modify temp != fix group");
    if (strcmp(temperature->style,"region") == 0 && iregion == -1 && 
	comm->me == 0)
      error->warning("Temperature for temp/rescale is style region");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::thermo_fields(int n, int *flags, char **keywords)
{
  if (n == 0) return 1;
  flags[0] = 3;
  strcpy(keywords[0],"EngTRscl");
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::thermo_compute(double *values)
{
  values[0] = energy;
  return 1;
}
