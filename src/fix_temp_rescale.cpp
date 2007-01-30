/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixTempRescale::FixTempRescale(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
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

  // create a new compute temp or temp/region style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[4];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  if (iregion == -1) {
    newarg[2] = "temp";
    modify->add_compute(3,newarg);
  } else {
    newarg[2] = "temp/region";
    newarg[3] = domain->regions[iregion]->id;
    modify->add_compute(4,newarg);
  }
  delete [] newarg;
  tflag = 1;
}

/* ---------------------------------------------------------------------- */

FixTempRescale::~FixTempRescale()
{
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::init()
{
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) error->all("Temp ID for fix temp/rescale does not exist");
  temperature = modify->compute[icompute];

  temperature->init();              // not yet called by Modify::init()
  efactor = (0.5 * force->boltz * temperature->dof);
  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixTempRescale::end_of_step()
{
  double t_current = temperature->compute_scalar();
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
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all("Could not find fix_modify temp ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all("Fix_modify temp ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning("Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixTempRescale::thermo(int n)
{
  if (n == 0) return energy;
  else return 0.0;
}
