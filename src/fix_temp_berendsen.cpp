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

#include "fix_temp_berendsen.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempBerendsen::FixTempBerendsen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  tstr(nullptr), id_temp(nullptr), tflag(0)
{
  if (narg != 6)
    error->all(FLERR,"Illegal fix {} command: expected 6 arguments but found {}", style, narg);

  // Berendsen thermostat should be applied every step

  restart_global = 1;
  dynamic_group_allow = 1;
  scalar_flag = 1;
  extscalar = 1;
  ecouple_flag = 1;
  nevery = 1;
  global_freq = nevery;

  tstr = nullptr;
  if (utils::strmatch(arg[3],"^v_")) {
    tstr = utils::strdup(arg[3]+2);
    tstyle = EQUAL;
  } else {
    t_start = utils::numeric(FLERR,arg[3],false,lmp);
    t_target = t_start;
    tstyle = CONSTANT;
  }

  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  t_period = utils::numeric(FLERR,arg[5],false,lmp);

  // error checks

  if (t_period <= 0.0)
    error->all(FLERR,"Fix temp/berendsen Tdamp period must be > 0.0");

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  tflag = 1;

  energy = 0;
}

/* ---------------------------------------------------------------------- */

FixTempBerendsen::~FixTempBerendsen()
{
  delete[] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete[] id_temp;
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
  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name {} for fix temp/berendsen does not exist", tstr);
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR,"Variable {} for fix temp/berendsen is invalid style", tstr);
  }

  temperature = modify->get_compute_by_id(id_temp);
  if (!temperature)
    error->all(FLERR,"Temperature compute ID {} for fix {} does not exist", id_temp, style);

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Cannot thermostat atoms in rigid bodies with fix {}", style);

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsen::end_of_step()
{
  double t_current = temperature->compute_scalar();
  double tdof = temperature->dof;

  // there is nothing to do, if there are no degrees of freedom

  if (tdof < 1) return;

  if (t_current == 0.0)
    error->all(FLERR, "Computed current temperature for fix temp/berendsen must not be 0.0");

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
      error->one(FLERR, "Fix temp/berendsen variable {} returned negative temperature",
                 input->variable->names[tvar]);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocities by lamda
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since lamda is multiplied by v

  double lamda = sqrt(1.0 + update->dt/t_period*(t_target/t_current - 1.0));
  double efactor = 0.5 * force->boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

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
  } else {
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
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR,"Could not find fix_modify temperature compute {}", id_temp);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature compute {} does not compute temperature", id_temp);
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group: {} vs {}",
                     group->names[igroup], group->names[temperature->igroup]);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsen::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempBerendsen::compute_scalar()
{
  return energy;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTempBerendsen::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = energy;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTempBerendsen::restart(char *buf)
{
  auto list = (double *) buf;

  energy = list[0];
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixTempBerendsen::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  }
  return nullptr;
}
