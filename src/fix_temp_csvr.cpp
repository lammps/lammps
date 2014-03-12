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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (ICTP, Italy)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_csvr.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempCSVR::FixTempCSVR(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix temp/csvr command");

  // CSVR thermostat should be applied every step

  nevery = 1;
  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;

  tstr = NULL;
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[3][2]);
    tstyle = EQUAL;
  } else {
    t_start = force->numeric(FLERR,arg[3]);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = force->numeric(FLERR,arg[4]);
  t_period = force->numeric(FLERR,arg[5]);
  int seed = force->inumeric(FLERR,arg[6]);

  // error checks

  if (t_period <= 0.0) error->all(FLERR,"Fix temp/csvr period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix temp/csvr random seed");

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

  random = new RanMars(lmp,seed + comm->me);

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempCSVR::~FixTempCSVR()
{
  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixTempCSVR::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempCSVR::init()
{
  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix temp/csvr does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR,"Variable for fix temp/csvr is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/csvr does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempCSVR::end_of_step()
{
  double t_current = temperature->compute_scalar();

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix temp/csvr variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // Langevin thermostat, implemented as decribed in
  // Bussi and Parrinello, Phys. Rev. E (2007).
  // it is a linear combination of old velocities and new,
  // randomly chosen, velocity, with proper coefficients

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  const double c1 = exp(-update->dt/t_period);

  if (atom->rmass_flag) { // per atom masses
    const double * const rmass = atom->rmass;
  
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          const double m = rmass[i];
          const double c2 = sqrt((1.0-c1*c1)*t_target/m);
          for (int j = 0; j < 3; ++j) {
            energy  += 0.5*m*v[i][j]*v[i][j];
            v[i][j] *= c1;
            v[i][j] += c2*random->gaussian();
            energy  -= 0.5*m*v[i][j]*v[i][j];
          }
        }
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          const double m = rmass[i];
          const double c2 = sqrt((1.0-c1*c1)*t_target/m);
          temperature->remove_bias(i,v[i]);
          for (int j = 0; j < 3; ++j) {
            energy  += 0.5*rmass[i]*v[i][j]*v[i][j];
            v[i][j] *= c1;
            v[i][j] += c2*random->gaussian();
            energy  -= 0.5*rmass[i]*v[i][j]*v[i][j];
          }
          temperature->restore_bias(i,v[i]);
        }
      }
    }
  } else { // per atom type masses

    const double * const mass = atom->mass;
    const int    * const type = atom->type;

    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          const double m = mass[type[i]];
          const double c2 = sqrt((1.0-c1*c1)*t_target/m);

          for (int j = 0; j < 3; ++j) {
            energy  += 0.5*m*v[i][j]*v[i][j];
            v[i][j] *= c1;
            v[i][j] += c2*random->gaussian();
            energy  -= 0.5*m*v[i][j]*v[i][j];
          }
        }
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          const double m = mass[type[i]];
          const double c2 = sqrt((1.0-c1*c1)*t_target/m);

          temperature->remove_bias(i,v[i]);
          for (int j = 0; j < 3; ++j) {
            energy  += 0.5*m*v[i][j]*v[i][j];
            v[i][j] *= c1;
            v[i][j] += c2*random->gaussian();
            energy  -= 0.5*m*v[i][j]*v[i][j];
          }
          temperature->restore_bias(i,v[i]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixTempCSVR::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempCSVR::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempCSVR::compute_scalar()
{
  return energy;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixTempCSVR::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  }
  return NULL;
}
