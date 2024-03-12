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

#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "math_const.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { NOCOUPLE = 0, XYZ, XY, YZ, XZ };

/* ---------------------------------------------------------------------- */

FixDeformPressure::FixDeformPressure(LAMMPS *lmp, int narg, char **arg) :
  FixDeform(lmp, narg, arg), id_temp(nullptr), id_press(nullptr), temperature(nullptr),
  pressure(nullptr)
{
  // set defaults

  set_extra = new SetExtra[7];
  memset(set_extra, 0, 7 * sizeof(SetExtra));
  memset(&set_box, 0, sizeof(Set));

  // parse only parameter/style arguments specific to this child class

  int index, iarg;
  std::size_t i = 0;
  while (i < leftover_iarg.size()) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "x") == 0 ||
        strcmp(arg[iarg], "y") == 0 ||
        strcmp(arg[iarg], "z") == 0) {

      if (strcmp(arg[iarg], "x") == 0) index = 0;
      else if (strcmp(arg[iarg], "y") == 0) index = 1;
      else if (strcmp(arg[iarg], "z") == 0) index = 2;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure", error);
      if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set_extra[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set_extra[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set_extra[index].pvar_flag = 1;
        }
        set_extra[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        i += 4;
      } else if (strcmp(arg[iarg + 1], "pressure/mean") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure/mean", error);
        set[index].style = PMEAN;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set_extra[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set_extra[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set_extra[index].pvar_flag = 1;
        }
        set_extra[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        i += 4;
      } else error->all(FLERR, "Illegal fix deform/pressure command argument: {}", arg[iarg + 1]);

    } else if (strcmp(arg[iarg], "xy") == 0 ||
               strcmp(arg[iarg], "xz") == 0 ||
               strcmp(arg[iarg], "yz") == 0) {

      if (strcmp(arg[iarg], "xy") == 0) index = 5;
      else if (strcmp(arg[iarg], "xz") == 0) index = 4;
      else if (strcmp(arg[iarg], "yz") == 0) index = 3;

      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure", error);
      if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set_extra[index].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set_extra[index].pstr = utils::strdup(&arg[iarg + 2][2]);
          set_extra[index].pvar_flag = 1;
        }
        set_extra[index].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        i += 4;
      } else error->all(FLERR, "Illegal fix deform/pressure command: {}", arg[iarg + 1]);

    } else if (strcmp(arg[iarg], "box") == 0) {
      if (strcmp(arg[iarg + 1], "volume") == 0) {
        set_box.style = VOLUME;
        i += 2;
      } else if (strcmp(arg[iarg + 1], "pressure") == 0) {
        if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure pressure", error);
        set_box.style = PRESSURE;
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2]) {
          set_extra[6].ptarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        } else {
          set_extra[6].pstr = utils::strdup(&arg[iarg + 2][2]);
          set_extra[6].pvar_flag = 1;
        }
        set_extra[6].pgain = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        i += 4;
      } else error->all(FLERR, "Illegal fix deform/pressure command argument: {}", arg[iarg + 1]);
    } else break;
  }

  // read options from end of input line
  // shift arguments before reading

  iarg = iarg_options_start;
  options(i, narg - iarg, &arg[iarg]);

  // repeat: setup dimflags used by other classes to check for volume-change conflicts

  for (int i = 0; i < 6; i++)
    if (set[i].style == NONE) dimflag[i] = 0;
    else dimflag[i] = 1;

  if (set_box.style != NONE) {
    dimflag[0] = 1;
    dimflag[1] = 1;
    dimflag[2] = 1;
  }

  if (dimflag[0]) box_change |= BOX_CHANGE_X;
  if (dimflag[1]) box_change |= BOX_CHANGE_Y;
  if (dimflag[2]) box_change |= BOX_CHANGE_Z;
  if (dimflag[3]) box_change |= BOX_CHANGE_YZ;
  if (dimflag[4]) box_change |= BOX_CHANGE_XZ;
  if (dimflag[5]) box_change |= BOX_CHANGE_XY;

  // repeat: no tensile deformation on shrink-wrapped dims
  // b/c shrink wrap will change box-length

  for (int i = 0; i < 3; i++)
    if (set_box.style && (domain->boundary[i][0] >= 2 || domain->boundary[i][1] >= 2))
      error->all(FLERR, "Cannot use fix deform/pressure on a shrink-wrapped boundary");

  // repeat: no tilt deformation on shrink-wrapped 2nd dim
  // b/c shrink wrap will change tilt factor in domain::reset_box()

  if (set[3].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR, "Cannot use fix deform/pressure tilt on a shrink-wrapped 2nd dim");
  if (set[4].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR, "Cannot use fix deform/pressure tilt on a shrink-wrapped 2nd dim");
  if (set[5].style && (domain->boundary[1][0] >= 2 || domain->boundary[1][1] >= 2))
    error->all(FLERR, "Cannot use fix deform/pressure tilt on a shrink-wrapped 2nd dim");

  // for VOLUME, setup links to other dims
  // fixed, dynamic1, dynamic2

  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;
    int other1 = (i + 1) % 3;
    int other2 = (i + 2) % 3;

    // Cannot use VOLUME option without at least one deformed dimension
    if (set[other1].style == NONE || set[other1].style == VOLUME)
      if (set[other2].style == NONE || set[other2].style == VOLUME)
        error->all(FLERR, "Fix {} volume setting is invalid", style);

    if (set[other1].style == NONE) {
      set[i].substyle = ONE_FROM_ONE;
      set[i].fixed = other1;
      set[i].dynamic1 = other2;
    } else if (set[other2].style == NONE) {
      set[i].substyle = ONE_FROM_ONE;
      set[i].fixed = other2;
      set[i].dynamic1 = other1;
    } else if (set[other1].style == VOLUME) {
      set[i].substyle = TWO_FROM_ONE;
      set[i].fixed = other1;
      set[i].dynamic1 = other2;
    } else if (set[other2].style == VOLUME) {
      set[i].substyle = TWO_FROM_ONE;
      set[i].fixed = other2;
      set[i].dynamic1 = other1;
    } else {
      set[i].substyle = ONE_FROM_TWO;
      set[i].dynamic1 = other1;
      set[i].dynamic2 = other2;
    }
  }

  // repeat: set varflag

  for (int i = 0; i < 7; i++)
    if (set_extra[i].pvar_flag) varflag = 1;

  // repeat: reneighboring only forced if flips can occur due to shape changes

  if ((!force_reneighbor) && flipflag && (set[3].style || set[4].style || set[5].style)) {
    force_reneighbor = 1;
    irregular = new Irregular(lmp);
  }

  // set initial values at time fix deform/pressure is issued

  set_box.vol_initial = domain->xprd * domain->yprd * domain->zprd;

  // populate coupled pressure controls

  if (pcouple != NOCOUPLE) {

    if (domain->dimension == 2)
      if (pcouple == XYZ || pcouple == XZ || pcouple == YZ)
        error->all(FLERR, "Cannot couple Z dimension in fix deform/pressure in 2D");

    int coupled_indices[3] = {0};
    int j = -1;

    if (pcouple == XYZ || pcouple == XY || pcouple == XZ)
      coupled_indices[0] = 1;
    if (pcouple == XYZ || pcouple == XY || pcouple == YZ)
      coupled_indices[1] = 1;
    if (pcouple == XYZ || pcouple == XZ || pcouple == YZ)
      coupled_indices[2] = 1;

    // Check coupled styles and find reference
    for (int i = 0; i < 3; i++) {
      if (coupled_indices[i]) {
        set_extra[i].coupled_flag = 1;
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
          set_extra[i].pgain = set_extra[j].pgain;
          if (set_extra[j].pvar_flag) {
            set_extra[i].pstr = set_extra[j].pstr;
            set_extra[i].pvar_flag = 1;
          } else {
            set_extra[i].ptarget = set_extra[j].ptarget;
          }
        } else {
          // Check for incompatibilities in style
          if (set[j].style != set[i].style && set[i].style != NONE)
            error->all(FLERR, "Cannot couple dimensions with different control options");
          if (set[j].style != PRESSURE) continue;

          // If pressure controlled, check for incompatibilities in parameters
          if (set_extra[i].pgain != set_extra[j].pgain || set_extra[i].pvar_flag != set_extra[j].pvar_flag ||
              set_extra[i].ptarget != set_extra[j].ptarget)
            error->all(FLERR, "Coupled dimensions must have identical gain parameters");

          if (set_extra[j].pvar_flag)
            if (strcmp(set_extra[i].pstr, set_extra[j].pstr) != 0)
              error->all(FLERR, "Coupled dimensions must have the same target pressure");
        }
      }
    }
  }

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
  for (int i = 0; i < 6; i++) {
    if (set[i].style == PRESSURE || set[i].style == PMEAN) {
      pressure_flag = 1;
      if (set_extra[i].pgain <= 0.0)
        error->all(FLERR, "Illegal fix deform/pressure gain constant, must be positive");
    }
    if (set_extra[i].coupled_flag) pressure_flag = 1;
  }
  if (set_box.style == PRESSURE) pressure_flag = 1;
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

  // check conflicts between x,y,z styles and box

  if (set_box.style)
    for (int i = 0; i < 3; i++)
      if (set[i].style == FINAL || set[i].style == DELTA || set[i].style == SCALE || set[i].style == PMEAN || set[i].style == VARIABLE)
        error->all(FLERR, "Cannot use fix deform/pressure box parameter with x, y, or z styles other than vel, erate, trate, pressure, and wiggle");

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
    temperature = modify->add_compute(fmt::format("{} all temp", id_temp));
    tflag = 1;

    // create a new compute pressure style
    // id = fix-ID + press, compute group = all
    // pass id_temp as 4th arg to pressure constructor

    id_press = utils::strdup(std::string(id) + "_press");
    pressure = modify->add_compute(fmt::format("{} all pressure {}", id_press, id_temp));
    pflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

FixDeformPressure::~FixDeformPressure()
{
  if (set_extra)
    for (int i = 0; i < 7; i++)
      delete[] set_extra[i].pstr;
  delete[] set_extra;

  delete[] set_box.hstr;
  delete[] set_box.hratestr;

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

  set_box.vol_start = domain->xprd * domain->yprd * domain->zprd;

  // check optional variables for PRESSURE or PMEAN style

  for (int i = 0; i < 7; i++) {
    if (!set_extra[i].pvar_flag) continue;
    set_extra[i].pvar = input->variable->find(set_extra[i].pstr);
    if (set_extra[i].pvar < 0)
      error->all(FLERR, "Variable name {} for fix deform/pressure does not exist", set_extra[i].pstr);
    if (!input->variable->equalstyle(set_extra[i].pvar))
      error->all(FLERR, "Variable {} for fix deform/pressure is invalid style", set_extra[i].pstr);
  }

  // Find pressure/temp computes if needed

  if (pressure_flag) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR, "Temperature ID {} for fix deform/pressure does not exist", id_temp);

    pressure = modify->get_compute_by_id(id_press);
    if (!pressure)
      error->all(FLERR, "Pressure ID {} for fix deform/pressure does not exist", id_press);
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

  if (strain_flag) FixDeform::apply_strain();

  // set new box size for pressure-based dims

  if (pressure_flag) {
    temperature->compute_vector();
    pressure->compute_vector();
    pressure->compute_scalar();
    for (int i = 0; i < 3; i++) {
      if (!set_extra[i].saved) {
        set_extra[i].saved = 1;
        set_extra[i].prior_rate = 0.0;
        set_extra[i].prior_pressure = pressure->vector[i];
      }
    }
    apply_pressure();
  }

  // set new box size for VOLUME dims that are linked to other dims
  // NOTE: still need to set h_rate for these dims

  if (volume_flag) apply_volume();

  // apply any final box scalings

  if (set_box.style) apply_box();

  // Save pressure/strain rate if required

  if (pressure_flag) {
    for (int i = 0; i < 3; i++) {
      set_extra[i].prior_pressure = pressure->vector[i];
      set_extra[i].prior_rate = ((set[i].hi_target - set[i].lo_target) /
                           (domain->boxhi[i] - domain->boxlo[i]) - 1.0)  / update->dt;
    }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);


  FixDeform::update_domain();

  // trigger virial computation, if needed, on next timestep

  if (pressure_flag)
    pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   apply pressure controls
------------------------------------------------------------------------- */

void FixDeformPressure::apply_pressure()
{
  // If variable pressure, calculate current target
  for (int i = 0; i < 6; i++)
    if (set[i].style == PRESSURE)
      if (set_extra[i].pvar_flag)
        set_extra[i].ptarget = input->variable->compute_equal(set_extra[i].pvar);

  // Find current (possibly coupled/hydrostatic) pressure for X, Y, Z
  double *tensor = pressure->vector;
  double scalar = pressure->scalar;
  double p_current[3] = {0.0, 0.0, 0.0};

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

    h_rate[i] = set_extra[i].pgain * (p_current[i] - set_extra[i].ptarget);

    if (normalize_pressure_flag) {
      if (set_extra[i].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);
      } else h_rate[i] /= fabs(set_extra[i].ptarget);
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

    h_rate[i] = L * set_extra[i].pgain * (pcurrent - set_extra[i].ptarget);
    if (normalize_pressure_flag) {
      if (set_extra[i].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);
      } else h_rate[i] /= fabs(set_extra[i].ptarget);
    }

    if (max_h_rate != 0)
      if (fabs(h_rate[i]) > max_h_rate)
        h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);

    set[i].tilt_target = tilt + update->dt * h_rate[i];
  }
}

/* ----------------------------------------------------------------------
   apply volume controls
------------------------------------------------------------------------- */

void FixDeformPressure::apply_volume()
{
  double e1, e2;
  int linked_pressure = 0;

  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;

    int dynamic1 = set[i].dynamic1;
    int dynamic2 = set[i].dynamic2;
    int fixed = set[i].fixed;
    double v0 = set[i].vol_start;
    double shift = 0.0;

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
        double e1i = set_extra[i].prior_rate;
        double e2i = set_extra[fixed].prior_rate;
        double L1i = domain->boxhi[i] - domain->boxlo[i];
        double L2i = domain->boxhi[fixed] - domain->boxlo[fixed];
        double L3i = domain->boxhi[dynamic1] - domain->boxlo[dynamic1];
        double L3 = (set[dynamic1].hi_target - set[dynamic1].lo_target);
        double Vi = L1i * L2i * L3i;
        double V = L3 * L1i * L2i;
        double e3 = (L3 / L3i - 1.0) / dt;
        double p1 = pressure->vector[i];
        double p2 = pressure->vector[fixed];
        double p1i = set_extra[i].prior_pressure;
        double p2i = set_extra[fixed].prior_pressure;
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
            if (max_h_rate != 0) {
              if ((fabs(e1) > max_h_rate) || (fabs(e2) > max_h_rate)) {
                if (fabs(e1) > fabs(e2))
                  adjust_linked_rates(e1, e2, e3, Vi, V);
                else
                  adjust_linked_rates(e2, e1, e3, Vi, V);
              }
            }
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
   apply box controls
------------------------------------------------------------------------- */

void FixDeformPressure::apply_box()
{
  int i;
  double scale, shift = 0.0;
  double v_rate;

  if (set_box.style == VOLUME) {
    double v0 = set_box.vol_start;
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

  } else if (set_box.style == PRESSURE) {

    // If variable pressure, calculate current target
    if (set_extra[6].pvar_flag)
      set_extra[6].ptarget = input->variable->compute_equal(set_extra[6].pvar);

    v_rate = set_extra[6].pgain * (pressure->scalar - set_extra[6].ptarget);

    if (normalize_pressure_flag) {
      if (set_extra[6].ptarget == 0) {
        if (max_h_rate == 0) {
          error->all(FLERR, "Cannot normalize error for zero pressure without defining a max rate");
        } else v_rate = max_h_rate * v_rate / fabs(v_rate);
      } else v_rate /= fabs(set_extra[6].ptarget);
    }

    if (max_h_rate != 0)
      if (fabs(v_rate) > max_h_rate)
        v_rate = max_h_rate * v_rate / fabs(v_rate);

    set_extra[6].cumulative_strain += update->dt * v_rate;
    scale = (1.0 + set_extra[6].cumulative_strain);
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

/* ----------------------------------------------------------------------
   write Set data to restart file
------------------------------------------------------------------------- */

void FixDeformPressure::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 9 * sizeof(double) + 7 * sizeof(Set) + 7 * sizeof(SetExtra);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(set, sizeof(Set), 6, fp);
    fwrite(&set_box, sizeof(Set), 1, fp);
    fwrite(set_extra, sizeof(SetExtra), 7, fp);
  }
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixDeformPressure::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;
  for (int i = 0; i < 6; i++)
    h_rate[i] = list[n++];
  for (int i = 0; i < 3; i++)
    h_ratelo[i] = list[n++];

  n = n * sizeof(double);
  int samestyle = 1;
  Set *set_restart = (Set *) &buf[n];
  for (int i = 0; i < 6; ++i) {
    // restore data from initial state
    set[i].lo_initial = set_restart[i].lo_initial;
    set[i].hi_initial = set_restart[i].hi_initial;
    set[i].vol_initial = set_restart[i].vol_initial;
    set[i].tilt_initial = set_restart[i].tilt_initial;
    // check if style settings are consistent (should do the whole set?)
    if (set[i].style != set_restart[i].style)
      samestyle = 0;
    if (set[i].substyle != set_restart[i].substyle)
      samestyle = 0;
  }
  n += 6 * sizeof(Set);

  // Only restore relevant box variables & check consistency
  Set set_box_restart;
  memcpy(&set_box_restart, (Set *) &buf[n], sizeof(Set));
  set_box.vol_initial = set_box_restart.vol_initial;
  if (set_box.style != set_box_restart.style)
    samestyle = 0;

  if (!samestyle)
    error->all(FLERR, "Fix deform/pressure settings not consistent with restart");

  n += sizeof(Set);
  SetExtra *set_extra_restart = (SetExtra *) &buf[n];
  for (int i = 0; i < 7; ++i) {
    set_extra[i].saved = set_extra_restart[i].saved;
    set_extra[i].prior_rate = set_extra_restart[i].prior_rate;
    set_extra[i].prior_pressure = set_extra_restart[i].prior_pressure;
    set_extra[i].cumulative_strain = set_extra_restart[i].cumulative_strain;
  }
}

/* ---------------------------------------------------------------------- */

void FixDeformPressure::options(int i, int narg, char **arg)
{
  pcouple = NOCOUPLE;
  max_h_rate = 0.0;
  vol_balance_flag = 0;
  normalize_pressure_flag = 0;

  // parse only options not handled by parent class

  int iarg;
  while (i < (int) leftover_iarg.size()) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "couple") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure couple", error);
      if (strcmp(arg[iarg + 1], "xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg + 1], "xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg + 1], "yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg + 1], "xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg + 1], "none") == 0) pcouple = NOCOUPLE;
      else error->all(FLERR, "Illegal fix deform/pressure couple command: {}", arg[iarg + 1]);
      i += 2;
    } else if (strcmp(arg[iarg], "max/rate") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure max/rate", error);
      max_h_rate = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (max_h_rate <= 0.0)
        error->all(FLERR, "Maximum strain rate must be a positive, non-zero value");
      i += 2;
    } else if (strcmp(arg[iarg], "normalize/pressure") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure normalize/pressure", error);
      normalize_pressure_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[iarg], "vol/balance/p") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix deform/pressure vol/balance/p", error);
      vol_balance_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 2;
    } else error->all(FLERR, "Illegal fix deform/pressure command: {}", arg[iarg]);
  }
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
