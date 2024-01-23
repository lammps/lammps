// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_deform_pressure.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "kspace.h"
#include "lattice.h"
#include "math_const.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE=0,FINAL,DELTA,SCALE,VEL,ERATE,TRATE,VOLUME,WIGGLE,VARIABLE,PRESSURE,PMEAN};
enum{ONE_FROM_ONE,ONE_FROM_TWO,TWO_FROM_ONE};
enum{NOCOUPLE=0,XYZ,XY,YZ,XZ};

/* ---------------------------------------------------------------------- */

FixDeformPressure::FixDeformPressure(LAMMPS *lmp, int narg, char **arg) : FixDeform(lmp, narg, arg),
id_temp(nullptr), id_press(nullptr)
{
  // populate coupled pressure controls

  if (pcouple != NOCOUPLE) {
    int coupled_indices[3] = {0};
    int j = -1;
    double couple_gain, coupled_pressure;
    char *couple_str;

    if (pcouple == XYZ || pcouple == XY || pcouple == XZ)
      coupled_indices[0] = 1;
    if (pcouple == XYZ || pcouple == XY || pcouple == YZ)
      coupled_indices[1] = 1;
    if (pcouple == XYZ || pcouple == XZ || pcouple == YZ)
      coupled_indices[2] = 1;

    // Check coupled styles and find reference
    for (int i = 0; i < 3; i++) {
      if (coupled_indices[i]) {
        set[i].coupled_flag = 1;
        if (set[i].style != PRESSURE && set[i].style != NONE)
          error->all(FLERR, "Cannot couple non-pressure-controlled dimensions");
        if (set[i].style == PRESSURE)
          j = i;
      }
    }

    if (j == -1)
      error->all(FLERR, "Must specify deformation style for at least one coupled dimension");

    // Copy or compare data for each coupled dimension
    for (int i = 0; i < 3; i++) {
      if (coupled_indices[i]) {
        // Copy coupling information if dimension style is undefined
        if (set[i].style == NONE) {
          set[i].style = PRESSURE;
          dimflag[i] = 1;
          set[i].pgain = set[j].pgain;
          if (set[j].pvar_flag) {
            set[i].pstr = set[j].pstr;
            set[i].pvar_flag = 1;
          } else {
            set[i].ptarget = set[j].ptarget;
          }
        } else {
          // Check for incompatibilities in style
          if (set[j].style != set[i].style && set[i].style != NONE)
            error->all(FLERR, "Cannot couple dimensions with different control options");
          if (set[j].style != PRESSURE) continue;

          // If pressure controlled, check for incompatibilities in parameters
          if (set[i].pgain != set[j].pgain || set[i].pvar_flag != set[j].pvar_flag ||
              set[i].ptarget != set[j].ptarget)
            error->all(FLERR, "Coupled dimensions must have identical gain parameters");

          if (set[j].pvar_flag)
            if (strcmp(set[i].pstr, set[j].pstr) != 0)
              error->all(FLERR, "Coupled dimensions must have the same target pressure");
        }
      }
    }
  }

  // repeat some checks in child class to catch changes to pcouple

  if (dimflag[0]) box_change |= BOX_CHANGE_X;
  if (dimflag[1]) box_change |= BOX_CHANGE_Y;
  if (dimflag[2]) box_change |= BOX_CHANGE_Z;

  // no tensile deformation on shrink-wrapped dims
  // b/c shrink wrap will change box-length

  for (int i = 0; i < 3; i++)
    if ((set[i].style || set[6].style) && (domain->boundary[i][0] >= 2 || domain->boundary[i][1] >= 2))
      error->all(FLERR, "Cannot use fix deform/pressure on a shrink-wrapped boundary");

  // if vol/balance/p used, must have 2 free dimensions

  if (vol_balance_flag) {
    for (int i = 0; i < 3; i++) {
      if (set[i].style != VOLUME) continue;
      if (set[i].substyle != TWO_FROM_ONE)
        error->all(FLERR, "Two dimensions must maintain constant volume to use the vol/balance/p option");
    }
  }

  // set strain_flag

  strain_flag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style != NONE && set[i].style != VOLUME &&
        set[i].style != PRESSURE && set[i].style != PMEAN)
      strain_flag = 1;

  // set pressure_flag

  pressure_flag = 0;
  for (int i = 0; i < 7; i++) {
    if (set[i].style == PRESSURE || set[i].style == PMEAN) pressure_flag = 1;
    if (set[i].coupled_flag) pressure_flag = 1;
  }
  if (vol_balance_flag) pressure_flag = 1;

  // check conflict between constant volume/pressure

  volume_flag = 0;
  for (int i = 0; i < 3; i++)
    if (set[i].style == VOLUME)
      volume_flag = 1;

  if (volume_flag)
    for (int i = 0; i < 6; i++)
      if (set[i].style == PMEAN)
        error->all(FLERR, "Cannot use fix deform/pressure to assign constant volume and pressure");

  // check conflicts between x,y,z styles and iso

  if (set[6].style)
    for (int i = 0; i < 3; i++)
      if (set[i].style == FINAL || set[i].style == DELTA || set[i].style == SCALE || set[i].style == PMEAN || set[i].style == VARIABLE)
        error->all(FLERR, "Cannot use fix deform/pressure iso parameter with x, y, or z styles other than vel, erate, trate, pressure, and wiggle");

  // check pressure used for max rate and normalize error flag

  if (!pressure_flag && max_h_rate != 0)
    error->all(FLERR, "Can only assign a maximum strain rate using pressure-controlled dimensions");

  if (!pressure_flag && normalize_pressure_flag)
    error->all(FLERR, "Can only normalize error using pressure-controlled dimensions");

  // Create pressure compute, if needed

  pflag = 0;
  tflag = 0;
  if (pressure_flag) {
    // create a new compute temp style
    // id = fix-ID + temp
    // compute group = all since pressure is always global (group all)
    //   and thus its KE/temperature contribution should use group all

    id_temp = utils::strdup(std::string(id) + "_temp");
    modify->add_compute(fmt::format("{} all temp",id_temp));
    tflag = 1;

    // create a new compute pressure style
    // id = fix-ID + press, compute group = all
    // pass id_temp as 4th arg to pressure constructor

    id_press = utils::strdup(std::string(id) + "_press");
    modify->add_compute(fmt::format("{} all pressure {}",id_press, id_temp));
    pflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

FixDeformPressure::~FixDeformPressure()
{
  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

void FixDeformPressure::init()
{
  FixDeform::init();

  // check optional variables for PRESSURE or PMEAN style

  for (int i = 0; i < 7; i++) {
    if (!set[i].pvar_flag) continue;
    set[i].pvar = input->variable->find(set[i].pstr);
    if (set[i].pvar < 0)
      error->all(FLERR, "Variable name {} for fix deform/pressure does not exist", set[i].pstr);
    if (!input->variable->equalstyle(set[i].pvar))
      error->all(FLERR, "Variable {} for fix deform/pressure is invalid style", set[i].pstr);
  }

  // Find pressure/temp computes if needed

  if (pressure_flag) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR, "Temperature ID for fix deform/pressure does not exist");
    temperature = modify->compute[icompute];

    icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all(FLERR, "Pressure ID for fix deform/pressure does not exist");
    pressure = modify->compute[icompute];
  }
}

/* ----------------------------------------------------------------------
   compute T,P if needed before integrator starts
------------------------------------------------------------------------- */

void FixDeformPressure::setup(int /*vflag*/)
{
  // trigger virial computation on next timestep
  if (pressure_flag) pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixDeformPressure::end_of_step()
{
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  // set new box size for strain-based dims

  if (strain_flag) FixDeform::set_strain();

  // set new box size for pressure-based dims

  if (pressure_flag) {
    temperature->compute_vector();
    pressure->compute_vector();
    pressure->compute_scalar();
    for (int i = 0; i < 3; i++) {
      if (!set[i].saved) {
        set[i].saved = 1;
        set[i].prior_rate = 0.0;
        set[i].prior_pressure = pressure->vector[i];
      }
    }
    set_pressure();
  }

  // set new box size for VOLUME dims that are linked to other dims
  // NOTE: still need to set h_rate for these dims

  if (volume_flag) set_volume();

  // apply any final isotropic scalings

  if (set[6].style) set_iso();

  // Save pressure/strain rate if required

  if (pressure_flag) {
    for (int i = 0; i < 3; i++) {
      set[i].prior_pressure = pressure->vector[i];
      set[i].prior_rate = ((set[i].hi_target - set[i].lo_target) /
                           (domain->boxhi[i] - domain->boxlo[i]) - 1.0)  / update->dt;
    }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);


  FixDeform::apply_deformation();

  // trigger virial computation, if needed, on next timestep

  if (pressure_flag)
    pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   set box size for pressure-based dimensions
------------------------------------------------------------------------- */

void FixDeformPressure::set_pressure()
{
  // If variable pressure, calculate current target
  for (int i = 0; i < 6; i++)
    if (set[i].style == PRESSURE)
      if (set[i].pvar_flag)
        set[i].ptarget = input->variable->compute_equal(set[i].pvar);

  // Find current (possibly coupled/hydrostatic) pressure for X, Y, Z
  double *tensor = pressure->vector;
  double scalar = pressure->scalar;
  double p_current[3];

  if (pcouple == XYZ) {
    double ave = THIRD * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else {
    if (set[0].style == PRESSURE) p_current[0] = tensor[0];
    else if (set[0].style == PMEAN) p_current[0] = scalar;

    if (set[1].style == PRESSURE) p_current[1] = tensor[1];
    else if (set[1].style == PMEAN) p_current[1] = scalar;

    if (set[2].style == PRESSURE) p_current[2] = tensor[2];
    else if (set[2].style == PMEAN) p_current[2] = scalar;
  }

  for (int i = 0; i < 3; i++) {
    if (set[i].style != PRESSURE && set[i].style != PMEAN) continue;

    h_rate[i] = set[i].pgain * (p_current[i] - set[i].ptarget);

    if (normalize_pressure_flag) {
      if (set[i].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);
      } else h_rate[i] /= fabs(set[i].ptarget);
    }

    if (max_h_rate != 0)
      if (fabs(h_rate[i]) > max_h_rate)
        h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);

    h_ratelo[i] = -0.5 * h_rate[i];

    double offset = 0.5 * (domain->boxhi[i] - domain->boxlo[i]) * (1.0 + update->dt * h_rate[i]);
    set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - offset;
    set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + offset;
  }

  for (int i = 3; i < 6; i++) {
    if (set[i].style != PRESSURE) continue;

    double L, tilt, pcurrent;
    if (i == 3) {
      L = domain->zprd;
      tilt = domain->yz;
      pcurrent = tensor[5];
    } else if (i == 4) {
      L = domain->zprd;
      tilt = domain->xz + update->dt;
      pcurrent = tensor[4];
    } else {
      L = domain->yprd;
      tilt = domain->xy;
      pcurrent = tensor[3];
    }

    h_rate[i] = L * set[i].pgain * (pcurrent - set[i].ptarget);
    if (normalize_pressure_flag) {
      if (set[i].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);
      } else h_rate[i] /= fabs(set[i].ptarget);
    }

    if (max_h_rate != 0)
      if (fabs(h_rate[i]) > max_h_rate)
        h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);

    set[i].tilt_target = tilt + update->dt * h_rate[i];
  }
}

/* ----------------------------------------------------------------------
   set box size for VOLUME dimensions
------------------------------------------------------------------------- */

void FixDeformPressure::set_volume()
{
  double e1, e2;
  int linked_pressure = 0;

  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;

    int dynamic1 = set[i].dynamic1;
    int dynamic2 = set[i].dynamic2;
    int fixed = set[i].fixed;
    double v0 = set[i].vol_start;
    double shift;

    if (set[i].substyle == ONE_FROM_ONE) {
      shift = 0.5 * (v0 / (set[dynamic1].hi_target - set[dynamic1].lo_target) /
             (set[fixed].hi_start-set[fixed].lo_start));
    } else if (set[i].substyle == ONE_FROM_TWO) {
      shift = 0.5 * (v0 / (set[dynamic1].hi_target - set[dynamic1].lo_target) /
             (set[dynamic2].hi_target - set[dynamic2].lo_target));
    } else if (set[i].substyle == TWO_FROM_ONE) {
      if (!vol_balance_flag) {
        shift = 0.5 * sqrt(v0 * (set[i].hi_start - set[i].lo_start) /
                   (set[dynamic1].hi_target - set[dynamic1].lo_target) /
                   (set[fixed].hi_start - set[fixed].lo_start));
      } else {
        double dt = update->dt;
        double e1i = set[i].prior_rate;
        double e2i = set[fixed].prior_rate;
        double L1i = domain->boxhi[i] - domain->boxlo[i];
        double L2i = domain->boxhi[fixed] - domain->boxlo[fixed];
        double L3i = domain->boxhi[dynamic1] - domain->boxlo[dynamic1];
        double L3 = (set[dynamic1].hi_target - set[dynamic1].lo_target);
        double Vi = L1i * L2i * L3i;
        double V = L3 * L1i * L2i;
        double e3 = (L3 / L3i - 1.0) / dt;
        double p1 = pressure->vector[i];
        double p2 = pressure->vector[fixed];
        double p1i = set[i].prior_pressure;
        double p2i = set[fixed].prior_pressure;
        double denominator;

        if (e3 == 0) {
          e1 = 0.0;
          e2 = 0.0;
          shift = 0.5 * L1i;
        } else if (e1i == 0 || e2i == 0 || (p2 == p2i && p1 == p1i)) {
          // If no prior strain or no change in pressure (initial step) just scale shift by relative box lengths
          shift = 0.5 * sqrt(v0 * L1i / L3 / L2i);
        } else {
          if (!linked_pressure) {
            // Calculate first strain rate by expanding stress to linear order, p1(t+dt) = p2(t+dt)
            // Calculate second strain rate to preserve volume
            denominator = p2 - p2i + e2i * ((p1 - p1i) / e1i);
            if (denominator != 0.0 && e1i != 0.0) {
              e1 = (((p2 - p2i) * (Vi - V) / (V * dt)) - e2i * (p1 - p2)) / denominator;
            } else {
              e1 = e2i;
            }
            e2 = (Vi - V * (1 + e1 * dt)) / (V * (1 + e1 * dt) * dt);

            // If strain rate exceeds limit in either dimension, cap it at the maximum compatible rate
            if (max_h_rate != 0)
              if (fabs(e1) > max_h_rate || fabs(e2) > max_h_rate)
                if (fabs(e1) > fabs(e2))
                  adjust_linked_rates(e1, e2, e3, Vi, V);
                else
                  adjust_linked_rates(e2, e1, e3, Vi, V);


            shift = 0.5 * L1i * (1.0 + e1 * dt);
            linked_pressure = 1;
          } else {
            // Already calculated value of e2
            shift = 0.5 * L1i * (1.0 + e2 * dt);
          }
        }
      }
    }

    h_rate[i] = (2.0 * shift / (domain->boxhi[i] - domain->boxlo[i]) - 1.0) / update->dt;
    h_ratelo[i] = -0.5 * h_rate[i];

    set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
    set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
  }
}


/* ----------------------------------------------------------------------
   Rescale volume preserving strain rates to enforce max rate
------------------------------------------------------------------------- */

void FixDeformPressure::adjust_linked_rates(double &e_larger, double &e_smaller, double e3, double Vi, double V)
{
  double dt = update->dt;
  double e_lim_positive = (Vi - V * (1 + max_h_rate * dt)) / (V * (1 + max_h_rate * dt) * dt);
  double e_lim_negative = (Vi - V * (1 - max_h_rate * dt)) / (V * (1 - max_h_rate * dt) * dt);
  if ((e_larger * e3) >= 0) {
    if (e_larger > 0.0) {
      // Same sign as primary strain rate, cap third dimension
      e_smaller = -max_h_rate;
      e_larger = e_lim_negative;
    } else {
      e_smaller = max_h_rate;
      e_larger = e_lim_positive;
    }
  } else {
    // Opposite sign, set to maxrate.
    if (e_larger > 0.0) {
      e_larger = max_h_rate;
      e_smaller = e_lim_positive;
    } else {
      e_larger = -max_h_rate;
      e_smaller = e_lim_negative;
    }
  }
}

/* ----------------------------------------------------------------------
   apply isotropic controls
------------------------------------------------------------------------- */

void FixDeformPressure::set_iso()
{
  int i;
  double scale, shift;
  double v_rate;

  if (set[6].style == VOLUME) {
    double v0 = set[6].vol_start;
    double v = 1.0;
    for (i = 0; i < 3; i++)
      v *= (set[i].hi_target - set[i].lo_target);

    scale = std::pow(v0 / v, THIRD);
    for (i = 0; i < 3; i++) {
      shift = 0.5 * (set[i].hi_target - set[i].lo_target) * scale;
      set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;

      // Recalculate h_rate
      h_rate[i] = (set[i].hi_target - set[i].lo_target) / (domain->boxhi[i] - domain->boxlo[i]) - 1.0;
      h_rate[i] /= update->dt;
      h_ratelo[i] = -0.5 * h_rate[i];
    }

  } else if (set[6].style == PRESSURE) {

    // If variable pressure, calculate current target
    if (set[6].pvar_flag)
      set[6].ptarget = input->variable->compute_equal(set[6].pvar);

    v_rate = set[6].pgain * (pressure->scalar- set[6].ptarget);

    if (normalize_pressure_flag) {
      if (set[6].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else v_rate = max_h_rate * v_rate / fabs(v_rate);
      } else v_rate /= fabs(set[6].ptarget);
    }

    if (max_h_rate != 0)
      if (fabs(v_rate) > max_h_rate)
        v_rate = max_h_rate * v_rate / fabs(v_rate);

    set[6].cumulative_strain += update->dt * v_rate;
    scale = (1.0 + set[6].cumulative_strain);
    for (i = 0; i < 3; i++) {
      shift = 0.5 * (set[i].hi_target - set[i].lo_target) * scale;
      set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;

      // Recalculate h_rate
      h_rate[i] = (set[i].hi_target - set[i].lo_target) / (domain->boxhi[i] - domain->boxlo[i]) - 1.0;
      h_rate[i] /= update->dt;
      h_ratelo[i] = -0.5 * h_rate[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDeformPressure::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal fix deform/pressure command");

  remapflag = Domain::X_REMAP;
  scaleflag = 1;
  flipflag = 1;

  pcouple = NOCOUPLE;
  dimension = domain->dimension;
  max_h_rate = 0.0;
  vol_balance_flag = 0;
  normalize_pressure_flag = 0;

  int index;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "x") == 0 ||
        strcmp(arg[iarg], "y") == 0 ||
        strcmp(arg[iarg], "z") == 0) {

      if (strcmp(arg[iarg], "x") == 0) index = 0;
      else if (strcmp(arg[iarg], "y") == 0) index = 1;
      else if (strcmp(arg[iarg], "z") == 0) index = 2;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure", error);
      if (strcmp(arg[iarg + 1], "final") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure final", error);
        set[index].style = FINAL;
        set[index].flo = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].fhi = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "delta") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure delta", error);
        set[index].style = DELTA;
        set[index].dlo = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].dhi = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "scale") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure scale", error);
        set[index].style = SCALE;
        set[index].scale = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "vel") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "erate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "trate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "volume") == 0) {
        set[index].style = VOLUME;
        iarg += 2;
      } else if (strcmp(arg[iarg + 1], "wiggle") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].tperiod = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure wiggle period, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "variable") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2])
          error->all(FLERR, "Illegal fix deform/pressure variable name {}", arg[iarg + 2]);
        if (strstr(arg[iarg + 3], "v_") != arg[iarg + 3])
          error->all(FLERR, "Illegal fix deform/pressure variable name {}", arg[iarg + 3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg + 2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg + 3][2]);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure pressure gain, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "pressure/mean") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure/mean", error);
        set[index].style = PMEAN;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure pressure gain, must be positive");
        iarg += 4;
      } else error->all(FLERR, "Illegal fix deform/pressure command argument: {}", arg[iarg + 1]);

    } else if (strcmp(arg[iarg], "xy") == 0 ||
               strcmp(arg[iarg], "xz") == 0 ||
               strcmp(arg[iarg], "yz") == 0) {

      if (triclinic == 0)
        error->all(FLERR, "fix deform/pressure tilt factors require triclinic box");
      if (strcmp(arg[iarg], "xy") == 0) index = 5;
      else if (strcmp(arg[iarg], "xz") == 0) index = 4;
      else if (strcmp(arg[iarg], "yz") == 0) index = 3;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure", error);
      if (strcmp(arg[iarg + 1], "final") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure final", error);
        set[index].style = FINAL;
        set[index].ftilt = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "delta") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure delta", error);
        set[index].style = DELTA;
        set[index].dtilt = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "vel") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "erate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "trate") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "wiggle") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        set[index].tperiod = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure wiggle period, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "variable") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2])
          error->all(FLERR, "Illegal fix deform/pressure variable name {}", arg[iarg + 2]);
        if (strstr(arg[iarg + 3], "v_") != arg[iarg + 3])
          error->all(FLERR, "Illegal fix deform/pressure variable name {}", arg[iarg + 3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg + 2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg + 3][2]);
        iarg += 4;
      } else if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure pressure gain, must be positive");
        iarg += 4;
      } else error->all(FLERR, "Illegal fix deform/pressure command: {}", arg[iarg + 1]);
    } else if (strcmp(arg[iarg], "iso") == 0) {
      index = 6;
      if (strcmp(arg[iarg + 1], "volume") == 0) {
        set[index].style = VOLUME;
        iarg += 2;
      } else if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR, "Illegal fix deform/pressure pressure gain, must be positive");
        iarg += 4;
      }
    } else break;
  }

  while (iarg < narg) {
    if (strcmp(arg[iarg], "remap") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure remap", error);
      if (strcmp(arg[iarg + 1], "x") == 0) remapflag = Domain::X_REMAP;
      else if (strcmp(arg[iarg + 1], "v") == 0) remapflag = Domain::V_REMAP;
      else if (strcmp(arg[iarg + 1], "none") == 0) remapflag = Domain::NO_REMAP;
      else error->all(FLERR, "Illegal fix deform/pressure remap command: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure units", error);
      if (strcmp(arg[iarg + 1], "box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0) scaleflag = 1;
      else error->all(FLERR, "Illegal fix deform/pressure units command: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "flip") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure flip", error);
      flipflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "couple") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure couple", error);
      if (strcmp(arg[iarg + 1], "xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg + 1], "xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg + 1], "yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg + 1], "xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg + 1], "none") == 0) pcouple = NOCOUPLE;
      else error->all(FLERR, "Illegal fix fix deform/pressure couple command: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "max/rate") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure max/rate", error);
      max_h_rate = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (max_h_rate <= 0.0)
        error->all(FLERR, "Maximum strain rate must be a positive, non-zero value");
      iarg += 2;
    } else if (strcmp(arg[iarg], "normalize/pressure") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure normalize/pressure", error);
      normalize_pressure_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "vol/balance/p") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure vol/balance/p", error);
      vol_balance_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else error->all(FLERR, "Illegal fix deform/pressure command: {}", arg[iarg]);
  }

  if (dimension == 2)
    if (pcouple == XYZ || pcouple == XZ || pcouple == YZ)
      error->all(FLERR, "Cannot couple Z dimension in fix deform/pressure in 2D");
}

/* ---------------------------------------------------------------------- */

int FixDeformPressure::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    temperature = modify->get_compute_by_id(arg[1]);
    if (!temperature)
      error->all(FLERR, "Could not find fix_modify temperature compute ID: ", arg[1]);

    if (temperature->tempflag == 0)
      error->all(FLERR, "Fix_modify temperature compute {} does not compute temperature", arg[1]);
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR, "Temperature compute {} for fix {} is not for group all: {}",
                     arg[1], style, group->names[temperature->igroup]);

    // reset id_temp of pressure to new temperature ID

    auto icompute = modify->get_compute_by_id(id_press);
    if (!icompute)
      error->all(FLERR, "Pressure compute ID {} for fix {} does not exist", id_press, style);
    icompute->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0], "press") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete[] id_press;
    id_press = utils::strdup(arg[1]);

    pressure = modify->get_compute_by_id(arg[1]);
    if (!pressure) error->all(FLERR, "Could not find fix_modify pressure compute ID: {}", arg[1]);
    if (pressure->pressflag == 0)
      error->all(FLERR, "Fix_modify pressure compute {} does not compute pressure", arg[1]);
    return 2;
  }

  return 0;
}
