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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_temp_csld.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempCSLD::FixTempCSLD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  vhold(nullptr), tstr(nullptr), id_temp(nullptr), random(nullptr)
{
  if (narg != 7) error->all(FLERR,"Illegal fix temp/csld command");

  // CSLD thermostat should be applied every step

  restart_global = 1;
  nevery = 1;
  scalar_flag = 1;
  ecouple_flag = 1;
  global_freq = nevery;
  dynamic_group_allow = 1;
  extscalar = 1;

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
  int seed = utils::inumeric(FLERR,arg[6],false,lmp);

  // error checks

  if (t_period <= 0.0) error->all(FLERR,"Illegal fix temp/csld command");
  if (seed <= 0) error->all(FLERR,"Illegal fix temp/csld  command");

  random = new RanMars(lmp,seed + comm->me);

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  tflag = 1;

  vhold = nullptr;
  nmax = -1;
  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempCSLD::~FixTempCSLD()
{
  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  delete random;
  memory->destroy(vhold);
  vhold = nullptr;
  nmax = -1;
}

/* ---------------------------------------------------------------------- */

int FixTempCSLD::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempCSLD::init()
{

  // we cannot handle constraints via rattle or shake correctly.

  int has_shake = 0;
  for (int i = 0; i < modify->nfix; i++)
    if ((strcmp(modify->fix[i]->style,"shake") == 0)
        || (strcmp(modify->fix[i]->style,"rattle") == 0)) ++has_shake;

  if (has_shake > 0)
    error->all(FLERR,"Fix temp/csld is not compatible with fix rattle or fix shake");

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix temp/csld does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else error->all(FLERR,"Variable for fix temp/csld is invalid style");
  }

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/csld does not exist");
  temperature = modify->compute[icompute];

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Cannot thermostat atoms in rigid bodies");

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempCSLD::end_of_step()
{

  // set current t_target
  // if variable temp, evaluate variable, wrap with clear/add

  double delta = update->ntimestep - update->beginstep;

  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix temp/csld variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  double t_current = temperature->compute_scalar();
  double ekin_old = t_current * 0.5 * temperature->dof * force->boltz;

  // there is nothing to do, if there are no degrees of freedom

  if (temperature->dof < 1) return;

  double * const * const v = atom->v;
  const int * const mask = atom->mask;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  // adjust holding space, if needed and copy existing velocities

  if (nmax < atom->nlocal) {
    nmax = atom->nlocal + 1;
    memory->destroy(vhold);
    memory->create(vhold,nmax,3,"csld:vhold");
  }

  // The CSLD thermostat is a linear combination of old and new velocities,
  // where the new ones are randomly chosen from a gaussian distribution.
  // see Bussi and Parrinello, Phys. Rev. E (2007).

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double m;
      if (atom->rmass_flag) m = atom->rmass[i];
      else m = atom->mass[type[i]];

      const double factor = 1.0/sqrt(m);
      const double vx = random->gaussian() * factor;
      vhold[i][0] = v[i][0];
      v[i][0] = vx;
      const double vy = random->gaussian() * factor;
      vhold[i][1] = v[i][1];
      v[i][1] = vy;
      const double vz = random->gaussian() * factor;
      vhold[i][2] = v[i][2];
      v[i][2] = vz;
    }
  }

  // mixing factors
  const double c1 = exp(-update->dt/t_period);
  const double c2 = sqrt((1.0-c1*c1)*t_target/temperature->compute_scalar());

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] = vhold[i][0]*c1 + v[i][0]*c2;
        v[i][1] = vhold[i][1]*c1 + v[i][1]*c2;
        v[i][2] = vhold[i][2]*c1 + v[i][2]*c2;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        temperature->remove_bias(i,vhold[i]);
        v[i][0] = vhold[i][0]*c1 + v[i][0]*c2;
        v[i][1] = vhold[i][1]*c1 + v[i][1]*c2;
        v[i][2] = vhold[i][2]*c1 + v[i][2]*c2;
        temperature->restore_bias(i,v[i]);
      }
    }
  }

  // tally the kinetic energy transferred between heat bath and system

  t_current = temperature->compute_scalar();
  energy +=  ekin_old - t_current * 0.5 * temperature->dof * force->boltz;
}

/* ---------------------------------------------------------------------- */

int FixTempCSLD::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

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

void FixTempCSLD::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempCSLD::compute_scalar()
{
  return energy;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTempCSLD::write_restart(FILE *fp)
{
  const int PRNGSIZE = 98+2+3;
  int nsize = PRNGSIZE*comm->nprocs+2; // pRNG state per proc + nprocs + energy
  double *list = nullptr;
  if (comm->me == 0) {
    list = new double[nsize];
    list[0] = energy;
    list[1] = comm->nprocs;
  }
  double state[PRNGSIZE];
  random->get_state(state);
  MPI_Gather(state,PRNGSIZE,MPI_DOUBLE,list+2,PRNGSIZE,MPI_DOUBLE,0,world);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
    delete[] list;
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTempCSLD::restart(char *buf)
{
  double *list = (double *) buf;

  energy = list[0];
  int nprocs = (int) list[1];
  if (nprocs != comm->nprocs) {
    if (comm->me == 0)
      error->warning(FLERR,"Different number of procs. Cannot restore RNG state.");
  } else random->set_state(list+2+comm->me*103);
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixTempCSLD::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  }
  return nullptr;
}
