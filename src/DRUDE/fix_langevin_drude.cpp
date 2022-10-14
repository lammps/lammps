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

#include "fix_langevin_drude.h"
#include "fix_drude.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixLangevinDrude::FixLangevinDrude(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix langevin/drude command");
  // TODO add option for tally

  // Langevin thermostat should be applied every step
  nevery = 1;
  global_freq = nevery;
  comm_reverse = 3;

  // core temperature
  tstr_core = nullptr;
  if (utils::strmatch(arg[3],"^v_")) {
    tstr_core = utils::strdup(arg[3]+2);
    tstyle_core = EQUAL;
  } else {
    t_start_core = utils::numeric(FLERR,arg[3],false,lmp);
    t_target_core = t_start_core;
    tstyle_core = CONSTANT;
  }
  t_period_core = utils::numeric(FLERR,arg[4],false,lmp);
  int seed_core = utils::inumeric(FLERR,arg[5],false,lmp);

  // drude temperature
  tstr_drude = nullptr;
  if (strstr(arg[7],"v_") == arg[6]) {
    tstr_drude = utils::strdup(arg[6]+2);
    tstyle_drude = EQUAL;
  } else {
    t_start_drude = utils::numeric(FLERR,arg[6],false,lmp);
    t_target_drude = t_start_drude;
    tstyle_drude = CONSTANT;
  }
  t_period_drude = utils::numeric(FLERR,arg[7],false,lmp);
  int seed_drude = utils::inumeric(FLERR,arg[8],false,lmp);

  // error checks
  if (t_period_core <= 0.0)
    error->all(FLERR,"Fix langevin/drude period must be > 0.0");
  if (seed_core  <= 0) error->all(FLERR,"Illegal langevin/drude seed");
  if (t_period_drude <= 0.0)
    error->all(FLERR,"Fix langevin/drude period must be > 0.0");
  if (seed_drude <= 0) error->all(FLERR,"Illegal langevin/drude seed");

  random_core  = new RanMars(lmp,seed_core);
  random_drude = new RanMars(lmp,seed_drude);

  int iarg = 9;
  zero = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"zero") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin/drude command");
      zero = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix langevin/drude command");
  }

  tflag = 0; // no external compute/temp is specified yet (for bias)
  energy = 0.;
  fix_drude = nullptr;
  temperature = nullptr;
  id_temp = nullptr;
}

/* ---------------------------------------------------------------------- */

FixLangevinDrude::~FixLangevinDrude()
{
  delete random_core;
  delete [] tstr_core;
  delete random_drude;
  delete [] tstr_drude;
}

/* ---------------------------------------------------------------------- */

int FixLangevinDrude::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::init()
{
  // check variable-style target core temperature
  if (tstr_core) {
    tvar_core = input->variable->find(tstr_core);
    if (tvar_core < 0)
      error->all(FLERR,"Variable name for fix langevin/drude does not exist");
    if (input->variable->equalstyle(tvar_core)) tstyle_core = EQUAL;
    else error->all(FLERR,"Variable for fix langevin/drude is invalid style");
  }

  // check variable-style target drude temperature
  if (tstr_drude) {
    tvar_drude = input->variable->find(tstr_drude);
    if (tvar_drude < 0)
      error->all(FLERR,"Variable name for fix langevin/drude does not exist");
    if (input->variable->equalstyle(tvar_drude)) tstyle_drude = EQUAL;
    else error->all(FLERR,"Variable for fix langevin/drude is invalid style");
  }

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"drude") == 0) break;
  if (ifix == modify->nfix) error->all(FLERR, "fix langevin/drude requires fix drude");
  fix_drude = dynamic_cast<FixDrude *>(modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::setup(int /*vflag*/)
{
  if (!utils::strmatch(update->integrate_style,"^verlet"))
    error->all(FLERR,"RESPA style not compatible with fix langevin/drude");
  if (!comm->ghost_velocity)
    error->all(FLERR,"fix langevin/drude requires ghost velocities. Use comm_modify vel yes");

  if (zero) {
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      int *drudetype = fix_drude->drudetype;
      int *type = atom->type;
      bigint ncore_loc = 0;
      for (int i=0; i<nlocal; i++)
          if (mask[i] & groupbit && drudetype[type[i]] != DRUDE_TYPE)
              ncore_loc++;
      MPI_Allreduce(&ncore_loc, &ncore, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  }
}

/* ---------------------------------------------------------------------- */

int FixLangevinDrude::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
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

void FixLangevinDrude::post_force(int /*vflag*/)
{
  // Thermalize by adding the langevin force if thermalize=true.
  // Each core-Drude pair is thermalized only once: where the core is local.

  double **v = atom->v, **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal, nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double ftm2v = force->ftm2v, mvv2e = force->mvv2e;
  double kb = force->boltz, dt = update->dt;

  int *drudetype = fix_drude->drudetype;
  tagint *drudeid = fix_drude->drudeid;
  double vdrude[3], vcore[3]; // velocities in reduced representation
  double fdrude[3], fcore[3]; // forces in reduced representation
  double Ccore, Cdrude, Gcore, Gdrude;
  double fcoresum[3], fcoreloc[3];
  int dim = domain->dimension;

  // Compute target core temperature
  if (tstyle_core == CONSTANT)
     t_target_core = t_start_core; // + delta * (t_stop-t_start_core);
  else {
      modify->clearstep_compute();
      t_target_core = input->variable->compute_equal(tvar_core);
      if (t_target_core < 0.0)
        error->one(FLERR, "Fix langevin/drude variable returned "
                          "negative core temperature");
      modify->addstep_compute(update->ntimestep + nevery);
  }

  // Compute target drude temperature
  if (tstyle_drude == CONSTANT)
      t_target_drude = t_start_drude; // + delta * (t_stop-t_start_core);
  else {
      modify->clearstep_compute();
      t_target_drude = input->variable->compute_equal(tvar_drude);
      if (t_target_drude < 0.0)
        error->one(FLERR, "Fix langevin/drude variable returned "
                          "negative drude temperature");
      modify->addstep_compute(update->ntimestep + nevery);
  }

  // Clear ghost forces
  // They have already been communicated if needed
  for (int i = nlocal; i < nall; i++) {
      for (int k = 0; k < dim; k++)
        f[i][k] = 0.;
  }
  if (zero) for (int k=0; k<dim; k++) fcoreloc[k] = 0.;

  // NB : the masses are the real masses, not the reduced ones
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { // only the cores need to be in the group
      if (drudetype[type[i]] == NOPOL_TYPE) { // Non-polarizable atom
        double mi;
        if (rmass)
          mi = rmass[i];
        else
          mi = mass[type[i]];
        Gcore  = mi / t_period_core  / ftm2v;
        Ccore  = sqrt(2.0 * Gcore  * kb * t_target_core  / dt / ftm2v / mvv2e);
        if (temperature) temperature->remove_bias(i, v[i]);
        for (int k = 0; k < dim; k++) {
            fcore[k] = Ccore  * random_core->gaussian()  - Gcore  * v[i][k];
            if (zero) fcoreloc[k] += fcore[k];
            f[i][k] += fcore[k];
        }
        if (temperature) temperature->restore_bias(i, v[i]);
      } else {
        if (drudetype[type[i]] == DRUDE_TYPE) continue; // do with the core

        int j = atom->map(drudeid[i]);
        double mi, mj, mtot, mu; // i is core, j is drude
        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          mi = mass[type[i]];
          mj = mass[type[j]];
        }
        mtot = mi + mj;
        mu = mi * mj / mtot;
        mi /= mtot;
        mj /= mtot;

        Gcore  = mtot / t_period_core  / ftm2v;
        Gdrude = mu   / t_period_drude / ftm2v;
        Ccore  = sqrt(2.0 * Gcore  * kb * t_target_core  / dt / ftm2v / mvv2e);
        Cdrude = sqrt(2.0 * Gdrude * kb * t_target_drude / dt / ftm2v / mvv2e);

        if (temperature) {
            temperature->remove_bias(i, v[i]);
            temperature->remove_bias(j, v[j]);
        }
        for (int k=0; k<dim; k++) {
          // TODO check whether a fix_modify temp can subtract a bias velocity
          vcore[k] = mi * v[i][k] + mj * v[j][k];
          vdrude[k] = v[j][k] - v[i][k];

          fcore[k]  = Ccore  * random_core->gaussian()  - Gcore  * vcore[k];
          fdrude[k] = Cdrude * random_drude->gaussian() - Gdrude * vdrude[k];

          if (zero) fcoreloc[k]  += fcore[k];

          f[i][k] += mi * fcore[k] - fdrude[k];
          f[j][k] += mj * fcore[k] + fdrude[k];

          // TODO tally energy if asked
        }
        if (temperature) {
            temperature->restore_bias(i, v[i]);
            temperature->restore_bias(j, v[j]);
        }
      }
    }
  }

  if (zero) { // Remove the drift
    MPI_Allreduce(fcoreloc, fcoresum, dim, MPI_DOUBLE, MPI_SUM, world);
    for (int k=0; k<dim; k++) fcoresum[k] /= ncore;
    for (int i=0; i<nlocal; i++) {
      if (mask[i] & groupbit) { // only the cores need to be in the group
        if (drudetype[type[i]] == NOPOL_TYPE) {
          for (int k=0; k<dim; k++) f[i][k] -= fcoresum[k];
        } else {
          if (drudetype[type[i]] == DRUDE_TYPE) continue; // do with the core
          int j = atom->map(drudeid[i]);
          double mi, mj, mtot; // i is core, j is drude
          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[type[i]];
            mj = mass[type[j]];
          }
          mtot = mi + mj;
          mi /= mtot;
          mj /= mtot;
          for (int k=0; k<dim; k++) {
            f[i][k] -= mi * fcoresum[k];
            f[j][k] -= mj * fcoresum[k];
          }
        }
      }
    }
  }

  // Reverse communication of the forces on ghost Drude particles
  comm->reverse_comm();
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::reset_target(double t_new)
{
  t_target_core = t_start_core = t_new;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixLangevinDrude::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"t_target_core") == 0) {
    return &t_target_core;
  } else if (strcmp(str,"t_target_drude") == 0) {
    return &t_target_drude;
  } else error->all(FLERR, "Illegal extract string in fix langevin/drude");
  return nullptr;
}

/* ---------------------------------------------------------------------- */
int FixLangevinDrude::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double ** f = atom->f;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixLangevinDrude::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double ** f = atom->f;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

