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

enum{STANDARD,REGION,PARTIAL};

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
  partial = 0;
  xflag = yflag = zflag = 1;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix temp/rescale command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) 
	error->all("Fix temp/rescale region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"partial") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix temp/rescale command");
      xflag = atoi(arg[iarg+1]);
      yflag = atoi(arg[iarg+2]);
      zflag = atoi(arg[iarg+3]);
      partial = 1;
      iarg += 4;
    } else error->all("Illegal fix temp/rescale command");
  }

  if (iregion == -1 && partial == 0) type = STANDARD;
  else if (iregion >= 0 && partial == 0) type = REGION;
  else if (iregion == -1 && partial) type = PARTIAL;
  else 
    error->all("Cannot use both region, partial options in fix temp/rescale");

  // create a new compute temp or temp/region or temp/partial
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[6];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  if (type == STANDARD) {
    newarg[2] = (char *) "temp";
    modify->add_compute(3,newarg);
  } else if (type == REGION) {
    newarg[2] = (char *) "temp/region";
    newarg[3] = domain->regions[iregion]->id;
    modify->add_compute(4,newarg);
  } else if (type == PARTIAL) {
    newarg[2] = (char *) "temp/partial";
    if (xflag) newarg[3] = (char *) "1";
    else newarg[3] = (char *) "0";
    if (yflag) newarg[4] = (char *) "1";
    else newarg[4] = (char *) "0";
    if (zflag) newarg[5] = (char *) "1";
    else newarg[5] = (char *) "0";
    modify->add_compute(6,newarg);
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

    if (type == STANDARD) {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
	}
      }

    } else if (type == REGION) {
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

    } else {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (xflag) v[i][0] *= factor;
	  if (yflag) v[i][1] *= factor;
	  if (zflag) v[i][2] *= factor;
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
