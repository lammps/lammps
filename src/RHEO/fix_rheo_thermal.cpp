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

#include "fix_rheo_thermal.h"
#include "fix_rheo.h"
#include "atom.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "update.h"
#include "force.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;
enum {NONE, CONSTANT, TYPE, ALUMINUM};

/* ---------------------------------------------------------------------- */

FixRHEOThermal::FixRHEOThermal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  Tc_style = NONE;
  cv_style = NONE;
  alpha_style = NONE;
  conductivity_style = NONE;
  dynamic_group_allow = 1;

  Tc_type = nullptr;
  kappa_type = nullptr;
  cv_type = nullptr;
  alpha_type = nullptr;

  int ntypes = atom->ntypes;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"conductivity") == 0) {
      // Conductivity arguments
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix command");
      if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix command");
        conductivity_style = CONSTANT;
        kappa = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (kappa < 0.0) error->all(FLERR,"Illegal fix command");
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+1+ntypes >= narg) error->all(FLERR,"Illegal fix command");
        conductivity_style = TYPE;
        memory->create(kappa_type,ntypes+1,"rheo_thermal:kappa_type");
        for (int i = 1; i <= ntypes; i++) {
          kappa_type[i] = utils::numeric(FLERR,arg[iarg+1+i],false,lmp);
          if (kappa_type[i] < 0.0) error->all(FLERR,"Illegal fix command");
        }
        iarg += 1+ntypes;
      } else {
        error->all(FLERR,"Illegal fix command");
      }
    } else if (strcmp(arg[iarg],"cv") == 0) {
      // Cv arguments
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix command");
      if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix command");
        cv_style = CONSTANT;
        cv = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (cv < 0.0) error->all(FLERR,"Illegal fix command");
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+1+ntypes >= narg) error->all(FLERR,"Illegal fix command");
        cv_style = TYPE;
        memory->create(cv_type,ntypes+1,"rheo_thermal:cv_type");
        for (int i = 1; i <= ntypes; i++) {
          cv_type[i] = utils::numeric(FLERR,arg[iarg+1+i],false,lmp);
          if (cv_type[i] < 0.0) error->all(FLERR,"Illegal fix command");
        }
        iarg += 1+ntypes;
      } else {
        error->all(FLERR,"Illegal fix command");
      }
    } else if (strcmp(arg[iarg],"alpha") == 0) {
      // Cv arguments
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix command");
      if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix command");
        alpha_style = CONSTANT;
        alpha = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+1+ntypes >= narg) error->all(FLERR,"Illegal fix command");
        alpha_style = TYPE;
        memory->create(alpha_type,ntypes+1,"rheo_thermal:alpha_type");
        for (int i = 1; i <= ntypes; i++) {
          alpha_type[i] = utils::numeric(FLERR,arg[iarg+1+i],false,lmp);
          if (alpha_type[i] < 0.0) error->all(FLERR,"Illegal fix command");
        }
        iarg += 1+ntypes;
      } else {
        error->all(FLERR,"Illegal fix command");
      }
    } else if (strcmp(arg[iarg],"Tfreeze") == 0) {
      // T freeze arguments
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix command");
      if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix command");
        Tc_style = CONSTANT;
        Tc = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (Tc < 0.0) error->all(FLERR,"Illegal fix command");
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+1+ntypes >= narg) error->all(FLERR,"Illegal fix command");
        Tc_style = TYPE;
        memory->create(Tc_type,ntypes+1,"rheo_thermal:Tc_type");
        for (int i = 1; i <= ntypes; i++) {
          Tc_type[i] = utils::numeric(FLERR,arg[iarg+1+i],false,lmp);
          if (Tc_type[i] < 0.0) error->all(FLERR,"Illegal fix command");
        }
        iarg += 1+ntypes;
      } else {
        error->all(FLERR,"Illegal fix command");
      }
    } else {
      error->all(FLERR,"Illegal fix command");
    }
    iarg += 1;
  }

  if (cv_style == NONE || conductivity_style == NONE)
    error->all(FLERR, "Must specify specific heat and conductivity styles\n");
}

/* ---------------------------------------------------------------------- */

FixRHEOThermal::~FixRHEOThermal()
{
  // If fix rheo is still defined, remove any set flags
  if (fix_rheo) {
    //fix_rheo->thermal_fix_defined = 0;
    //if(viscosity_style != NONE) fix_rheo->viscosity_fix_defined = 0;
  }

  memory->destroy(cv_type);
  memory->destroy(Tc_type);
  memory->destroy(kappa_type);
  memory->destroy(alpha_type);
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::init()
{
  int flag;
  int ifix = modify->find_fix_by_style("rheo");
  if (ifix == -1) error->all(FLERR, "Need to define fix rheo to use fix rheo/thermal");
  fix_rheo = ((FixRHEO *) modify->fix[ifix]);

  if (!fix_rheo->thermal_flag) error->all(FLERR, "Need to define thermal setting in fix rheo");

  if (fix_rheo->thermal_fix_defined)
    error->all(FLERR, "Cannot define two rheo thermal evolution fixes");
  fix_rheo->thermal_fix_defined = 1;

  int ifix2 = modify->find_fix_by_style("rheo/thermal");
  if (ifix > ifix2)
    error->all(FLERR, "Fix RHEO must be defined before fix RHEO/thermal");

  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::initial_integrate(int /*vflag*/)
{
  // update temperature from shifting
  int i;
  int *status = atom->status;
  double **gradt = compute_grad->gradt;
  double **vshift = compute_vshift->array_atom;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  if (shift_flag) {
    compute_vshift->correct_surfaces();
    for (i = 0; i < nlocal; i++) {

      if (!(status[i] & STATUS_SHIFT)) continue;

      if (mask[i] & groupbit) {
        for (a = 0; a < dim; a++) {
          temperature[i] += dtv * vshift[i][a] * gradt[i][a];
        }
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

void FixRHEOThermal::setup_pre_force(int /*vflag*/)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::post_integrate()
{
  //Post integrate teo ensure gradient up to date
  int *phase = atom->phase;
  double *temp = atom->temp;
  double *heat = atom->heat;
  double *rho = atom->rho;
  int *mask = atom->mask;

  double cvi, Tc, Tci, Ti, alphai;

  //Integrate temperature and check phase
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (phase[i] == FixRHEO::FLUID_NO_FORCE) continue;

      cvi = calc_cv(i);
      temp[i] += dtf*heat[i]/cvi;

      if (alpha_style != NONE && phase[i] <= FixRHEO::FLUID_MAX) {
        alphai = calc_alpha(i);
        rho[i] += dtf*heat[i]/cvi*alphai;
      }

      if (Tc_style != NONE) {
        Ti = temp[i];
        Tci = calc_Tc(i);

        if (Ti > Tci) {
          phase[i] = FixRHEO::FLUID;
        } else {
          if (phase[i] == FixRHEO::SOLID) {
            continue;
          } else {
            phase[i] = FixRHEO::FREEZING;
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::pre_force(int /*vflag*/)
{
  double *conductivity = atom->conductivity;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Calculate non-persistent quantities before pairstyles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      conductivity[i] = calc_kappa(i);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::end_of_step()
{
  double *temp = atom->temp;
  double *heat = atom->heat;
  int *phase = atom->phase;
  int *mask = atom->mask;
  double *rho = atom->rho;

  double cvi, alphai;

  //Integrate temperature and check phase
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (phase[i] == FixRHEO::FLUID_NO_FORCE) continue;

      cvi = calc_cv(i);
      temp[i] += dtf*heat[i]/cvi;

      if (alpha_style != NONE && phase[i] <= FixRHEO::FLUID_MAX) {
        alphai = calc_alpha(i);
        rho[i] += dtf*heat[i]/cvi*alphai;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_kappa(int i)
{
  if (conductivity_style == CONSTANT) {
    return kappa;
  } else if (conductivity_style == TYPE) {
    int itype = atom->type[i];
    return(kappa_type[itype]);
  } else {
    error->all(FLERR, "Invalid style");
  }
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_cv(int i)
{
  if (cv_style == CONSTANT) {
    return cv;
  } else if (cv_style == TYPE) {
    int itype = atom->type[i];
    return(cv_type[itype]);
  } else {
    error->all(FLERR, "Invalid style");
  }
}


/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_alpha(int i)
{
  if (alpha_style == CONSTANT) {
    return alpha;
  } else if (alpha_style == TYPE) {
    int itype = atom->type[i];
    return(alpha_type[itype]);
  } else {
    error->all(FLERR, "Invalid style");
  }
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_Tc(int i)
{
  if (Tc_style == CONSTANT) {
    return Tc;
  } else if (Tc_style == TYPE) {
    int itype = atom->type[i];
    return(Tc_type[itype]);
  } else {
    error->all(FLERR, "Invalid style");
  }
}