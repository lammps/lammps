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
#include "fix_temp_berendsen.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixTempBerendsen::FixTempBerendsen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix temp/berendsen command");

  // Berendsen thermostat should be applied every step

  nevery = 1;

  t_start = atof(arg[3]);
  t_stop = atof(arg[4]);
  t_period = atof(arg[5]);

  // error checks

  if (t_period <= 0.0) error->all("Fix temp/berendsen period must be > 0.0");

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;
}

/* ---------------------------------------------------------------------- */

FixTempBerendsen::~FixTempBerendsen()
{
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempBerendsen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsen::init()
{
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all("Temperature ID for fix temp/berendsen does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsen::end_of_step()
{
  double t_current = temperature->compute_scalar();
  if (t_current == 0.0)
    error->all("Computed temperature for fix temp/berendsen cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop-t_start);

  // rescale velocities by lamda

  double lamda = sqrt(1.0 + update->dt/t_period*(t_target/t_current - 1.0));

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	v[i][0] *= lamda;
	v[i][1] *= lamda;
	v[i][2] *= lamda;
      }
    }
  } else if (which == BIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	temperature->remove_bias(i,v[i]);
	v[i][0] *= lamda;
	v[i][1] *= lamda;
	v[i][2] *= lamda;
	temperature->restore_bias(i,v[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixTempBerendsen::modify_param(int narg, char **arg)
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
    if (icompute < 0) error->all("Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all("Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning("Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsen::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}
