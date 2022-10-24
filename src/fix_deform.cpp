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
   Contributing author: Pieter in 't Veld (SNL), Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_deform.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
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

FixDeform::FixDeform(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
irregular(nullptr), set(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix deform command");

  no_change_box = 1;
  restart_global = 1;
  pre_exchange_migrate = 1;
  pcouple = NOCOUPLE;
  dimension = domain->dimension;
  max_h_rate = 0.0;
  vol_balance_flag = 0;
  normalize_pressure_flag = 0;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix deform command");

  // set defaults

  set = new Set[6];
  memset(set,0,6*sizeof(Set));

  // parse arguments

  triclinic = domain->triclinic;

  int index;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0 ||
        strcmp(arg[iarg],"y") == 0 ||
        strcmp(arg[iarg],"z") == 0) {

      if (strcmp(arg[iarg],"x") == 0) index = 0;
      else if (strcmp(arg[iarg],"y") == 0) index = 1;
      else if (strcmp(arg[iarg],"z") == 0) index = 2;

      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform", error);
      if (strcmp(arg[iarg+1],"final") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform final", error);
        set[index].style = FINAL;
        set[index].flo = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        set[index].fhi = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform delta", error);
        set[index].style = DELTA;
        set[index].dlo = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        set[index].dhi = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"scale") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform scale", error);
        set[index].style = SCALE;
        set[index].scale = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"vel") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"erate") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"trate") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"volume") == 0) {
        set[index].style = VOLUME;
        iarg += 2;
      } else if (strcmp(arg[iarg+1],"wiggle") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        set[index].tperiod = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR,"Illegal fix deform wiggle period, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"variable") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2])
          error->all(FLERR,"Illegal fix deform variable name {}", arg[iarg+2]);
        if (strstr(arg[iarg+3],"v_") != arg[iarg+3])
          error->all(FLERR,"Illegal fix deform variable name {}", arg[iarg+3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg+2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg+3][2]);
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"pressure") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2]) {
          set[index].ptarget = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg+2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR,"Illegal fix deform pressure gain, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"pressure/mean") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform pressure/mean", error);
        set[index].style = PMEAN;
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2]) {
          set[index].ptarget = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg+2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR,"Illegal fix deform pressure gain, must be positive");
        iarg += 4;
      } else error->all(FLERR,"Illegal fix deform command argument: {}", arg[iarg+1]);

    } else if (strcmp(arg[iarg],"xy") == 0 ||
               strcmp(arg[iarg],"xz") == 0 ||
               strcmp(arg[iarg],"yz") == 0) {

      if (triclinic == 0)
        error->all(FLERR,"Fix deform tilt factors require triclinic box");
      if (strcmp(arg[iarg],"xy") == 0) index = 5;
      else if (strcmp(arg[iarg],"xz") == 0) index = 4;
      else if (strcmp(arg[iarg],"yz") == 0) index = 3;

      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform", error);
      if (strcmp(arg[iarg+1],"final") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform final", error);
        set[index].style = FINAL;
        set[index].ftilt = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform delta", error);
        set[index].style = DELTA;
        set[index].dtilt = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"vel") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform vel", error);
        set[index].style = VEL;
        set[index].vel = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"erate") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform erate", error);
        set[index].style = ERATE;
        set[index].rate = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"trate") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix deform trate", error);
        set[index].style = TRATE;
        set[index].rate = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"wiggle") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform wiggle", error);
        set[index].style = WIGGLE;
        set[index].amplitude = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        set[index].tperiod = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (set[index].tperiod <= 0.0)
          error->all(FLERR,"Illegal fix deform wiggle period, must be positive");
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"variable") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform variable", error);
        set[index].style = VARIABLE;
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2])
          error->all(FLERR,"Illegal fix deform variable name {}", arg[iarg+2]);
        if (strstr(arg[iarg+3],"v_") != arg[iarg+3])
          error->all(FLERR,"Illegal fix deform variable name {}", arg[iarg+3]);
        delete[] set[index].hstr;
        delete[] set[index].hratestr;
        set[index].hstr = utils::strdup(&arg[iarg+2][2]);
        set[index].hratestr = utils::strdup(&arg[iarg+3][2]);
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"pressure") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix deform pressure", error);
        set[index].style = PRESSURE;
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2]) {
          set[index].ptarget = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        } else {
          set[index].pstr = utils::strdup(&arg[iarg+2][2]);
          set[index].pvar_flag = 1;
        }
        set[index].pgain = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (set[index].pgain <= 0.0)
          error->all(FLERR,"Illegal fix deform pressure gain, must be positive");
        iarg += 4;
      } else error->all(FLERR,"Illegal fix deform command: {}", arg[iarg+1]);
    } else break;
  }

  // read options from end of input line
  // no x remap effectively moves atoms within box, so set restart_pbc

  options(narg-iarg,&arg[iarg]);
  if (remapflag != Domain::X_REMAP) restart_pbc = 1;

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

  // setup dimflags used by other classes to check for volume-change conflicts

  for (int i = 0; i < 6; i++)
    if (set[i].style == NONE) dimflag[i] = 0;
    else dimflag[i] = 1;

  if (dimflag[0]) box_change |= BOX_CHANGE_X;
  if (dimflag[1]) box_change |= BOX_CHANGE_Y;
  if (dimflag[2]) box_change |= BOX_CHANGE_Z;
  if (dimflag[3]) box_change |= BOX_CHANGE_YZ;
  if (dimflag[4]) box_change |= BOX_CHANGE_XZ;
  if (dimflag[5]) box_change |= BOX_CHANGE_XY;

  // no tensile deformation on shrink-wrapped dims
  // b/c shrink wrap will change box-length

  for (int i = 0; i < 3; i++)
    if (set[i].style && (domain->boundary[i][0] >= 2 || domain->boundary[i][1] >= 2))
      error->all(FLERR,"Cannot use fix deform on a shrink-wrapped boundary");

  // no tilt deformation on shrink-wrapped 2nd dim
  // b/c shrink wrap will change tilt factor in domain::reset_box()

  if (set[3].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR,"Cannot use fix deform tilt on a shrink-wrapped 2nd dim");
  if (set[4].style && (domain->boundary[2][0] >= 2 || domain->boundary[2][1] >= 2))
    error->all(FLERR,"Cannot use fix deform tilt on a shrink-wrapped 2nd dim");
  if (set[5].style && (domain->boundary[1][0] >= 2 || domain->boundary[1][1] >= 2))
    error->all(FLERR,"Cannot use fix deform tilt on a shrink-wrapped 2nd dim");

  // apply scaling to FINAL,DELTA,VEL,WIGGLE since they have dist/vel units

  int flag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == VEL || set[i].style == WIGGLE) flag = 1;

  double xscale,yscale,zscale;
  if (flag && scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // for 3,4,5: scaling is in 1st dimension, e.g. x for xz

  double map[6];
  map[0] = xscale; map[1] = yscale; map[2] = zscale;
  map[3] = yscale; map[4] = xscale; map[5] = xscale;

  for (int i = 0; i < 3; i++) {
    if (set[i].style == FINAL) {
      set[i].flo *= map[i];
      set[i].fhi *= map[i];
    } else if (set[i].style == DELTA) {
      set[i].dlo *= map[i];
      set[i].dhi *= map[i];
    } else if (set[i].style == VEL) {
      set[i].vel *= map[i];
    } else if (set[i].style == WIGGLE) {
      set[i].amplitude *= map[i];
    }
  }

  for (int i = 3; i < 6; i++) {
    if (set[i].style == FINAL) set[i].ftilt *= map[i];
    else if (set[i].style == DELTA) set[i].dtilt *= map[i];
    else if (set[i].style == VEL) set[i].vel *= map[i];
    else if (set[i].style == WIGGLE) set[i].amplitude *= map[i];
  }

  // for VOLUME, setup links to other dims
  // fixed, dynamic1, dynamic2

  volume_flag = 0;
  for (int i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;
    volume_flag = 1;
    int other1 = (i+1) % 3;
    int other2 = (i+2) % 3;

    // Cannot use VOLUME option without at least one deformed dimension
    if (set[other1].style == NONE || set[other1].style == VOLUME)
      if (set[other2].style == NONE || set[other2].style == VOLUME)
        error->all(FLERR,"Fix deform volume setting is invalid");

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

    if (vol_balance_flag && set[i].substyle != TWO_FROM_ONE)
      error->all(FLERR, "Two dimensions must maintain constant volume to use the vol/balance/p option");
  }

  // set strain_flag

  strain_flag = 0;
  for (int i = 0; i < 6; i++)
    if (set[i].style != NONE && set[i].style != VOLUME &&
        set[i].style != PRESSURE && set[i].style != PMEAN)
      strain_flag = 1;

  // set varflag

  varflag = 0;
  for (int i = 0; i < 6; i++) {
    if (set[i].style == VARIABLE) varflag = 1;
    if (set[i].pvar_flag) varflag = 1;
  }

  // set pressure_flag

  pressure_flag = 0;
  for (int i = 0; i < 6; i++) {
    if (set[i].style == PRESSURE || set[i].style == PMEAN) pressure_flag = 1;
    if (set[i].coupled_flag) pressure_flag = 1;
  }
  if (vol_balance_flag) pressure_flag = 1;

  // check conflict between constant volume/pressure

  if (volume_flag)
    for (int i = 0; i < 6; i++)
      if (set[i].style == PMEAN)
        error->all(FLERR, "Cannot use fix deform to assign constant volume and pressure");

  // check pressure used for max rate and normalize error flag

  if (!pressure_flag && max_h_rate != 0)
    error->all(FLERR, "Can only assign a maximum strain rate using pressure-controlled dimensions");

  if (!pressure_flag && normalize_pressure_flag)
    error->all(FLERR, "Can only normalize error using pressure-controlled dimensions");

  // set initial values at time fix deform is issued

  for (int i = 0; i < 3; i++) {
    set[i].lo_initial = domain->boxlo[i];
    set[i].hi_initial = domain->boxhi[i];
    set[i].vol_initial = domain->xprd * domain->yprd * domain->zprd;
  }
  set[3].tilt_initial = domain->yz;
  set[4].tilt_initial = domain->xz;
  set[5].tilt_initial = domain->xy;

  // reneighboring only forced if flips can occur due to shape changes

  if (flipflag && (set[3].style || set[4].style || set[5].style))
    force_reneighbor = 1;
  next_reneighbor = -1;

  flip = 0;

  if (force_reneighbor) irregular = new Irregular(lmp);
  else irregular = nullptr;

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

FixDeform::~FixDeform()
{
  if (set) {
    for (int i = 0; i < 6; i++) {
      delete[] set[i].hstr;
      delete[] set[i].hratestr;
      delete[] set[i].pstr;
    }
  }
  delete[] set;

  delete irregular;

  // reset domain's h_rate = 0.0, since this fix may have made it non-zero

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixDeform::setmask()
{
  int mask = 0;
  if (force_reneighbor) mask |= PRE_EXCHANGE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeform::init()
{
  // error if more than one fix deform
  // domain, fix nvt/sllod, compute temp/deform only work on single h_rate

  if (modify->get_fix_by_style("deform").size() > 1)
    error->all(FLERR,"More than one fix deform");

  // Kspace setting

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // elapsed time for entire simulation, including multiple runs if defined

  double delt = (update->endstep - update->beginstep) * update->dt;

  // check variables for VARIABLE style

  for (int i = 0; i < 6; i++) {
    if (set[i].style != VARIABLE) continue;
    set[i].hvar = input->variable->find(set[i].hstr);
    if (set[i].hvar < 0)
      error->all(FLERR,"Variable name {} for fix deform does not exist", set[i].hstr);
    if (!input->variable->equalstyle(set[i].hvar))
      error->all(FLERR,"Variable {} for fix deform is invalid style", set[i].hstr);
    set[i].hratevar = input->variable->find(set[i].hratestr);
    if (set[i].hratevar < 0)
      error->all(FLERR,"Variable name {} for fix deform does not exist", set[i].hratestr);
    if (!input->variable->equalstyle(set[i].hratevar))
      error->all(FLERR,"Variable {} for fix deform is invalid style", set[i].hratestr);
  }

  // check optional variables for PRESSURE or PMEAN style

  for (int i = 0; i < 6; i++) {
    if (!set[i].pvar_flag) continue;
    set[i].pvar = input->variable->find(set[i].pstr);
    if (set[i].pvar < 0)
      error->all(FLERR,"Variable name {} for fix deform does not exist", set[i].pstr);
    if (!input->variable->equalstyle(set[i].pvar))
      error->all(FLERR,"Variable {} for fix deform is invalid style", set[i].pstr);
  }

  // set start/stop values for box size and shape
  // if single run, start is current values
  // if multiple runs enabled via run start/stop settings,
  //   start is value when fix deform was issued
  // if VARIABLE or NONE, no need to set

  for (int i = 0; i < 3; i++) {
    if (update->firststep == update->beginstep) {
      set[i].lo_start = domain->boxlo[i];
      set[i].hi_start = domain->boxhi[i];
      set[i].vol_start = domain->xprd * domain->yprd * domain->zprd;
    } else {
      set[i].lo_start = set[i].lo_initial;
      set[i].hi_start = set[i].hi_initial;
      set[i].vol_start = set[i].vol_initial;
    }

    if (set[i].style == FINAL) {
      set[i].lo_stop = set[i].flo;
      set[i].hi_stop = set[i].fhi;
    } else if (set[i].style == DELTA) {
      set[i].lo_stop = set[i].lo_start + set[i].dlo;
      set[i].hi_stop = set[i].hi_start + set[i].dhi;
    } else if (set[i].style == SCALE) {
      double shift = 0.5 * set[i].scale * (set[i].hi_start - set[i].lo_start);
      set[i].lo_stop = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_stop = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
    } else if (set[i].style == VEL) {
      set[i].lo_stop = set[i].lo_start - 0.5 * delt * set[i].vel;
      set[i].hi_stop = set[i].hi_start + 0.5 * delt * set[i].vel;
    } else if (set[i].style == ERATE) {
      double shift = 0.5 * delt * set[i].rate * (set[i].hi_start - set[i].lo_start);
      set[i].lo_stop = set[i].lo_start - shift;
      set[i].hi_stop = set[i].hi_start + shift;
      if (set[i].hi_stop <= set[i].lo_stop)
        error->all(FLERR,"Final box dimension due to fix deform is < 0.0");
    } else if (set[i].style == TRATE) {
      double shift = 0.5 * ((set[i].hi_start - set[i].lo_start) * exp(set[i].rate * delt));
      set[i].lo_stop = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_stop = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
    } else if (set[i].style == WIGGLE) {
      double shift = 0.5 * set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].lo_stop = set[i].lo_start - shift;
      set[i].hi_stop = set[i].hi_start + shift;
    }
  }

  for (int i = 3; i < 6; i++) {
    if (update->firststep == update->beginstep) {
      if (i == 5) set[i].tilt_start = domain->xy;
      else if (i == 4) set[i].tilt_start = domain->xz;
      else if (i == 3) set[i].tilt_start = domain->yz;
    } else set[i].tilt_start = set[i].tilt_initial;

    if (set[i].style == FINAL) {
      set[i].tilt_stop = set[i].ftilt;
    } else if (set[i].style == DELTA) {
      set[i].tilt_stop = set[i].tilt_start + set[i].dtilt;
    } else if (set[i].style == VEL) {
      set[i].tilt_stop = set[i].tilt_start + delt * set[i].vel;
    } else if (set[i].style == ERATE) {
      if (i == 3) set[i].tilt_stop = set[i].tilt_start +
                    delt * set[i].rate * (set[2].hi_start - set[2].lo_start);
      if (i == 4) set[i].tilt_stop = set[i].tilt_start +
                    delt * set[i].rate * (set[2].hi_start - set[2].lo_start);
      if (i == 5) set[i].tilt_stop = set[i].tilt_start +
                    delt * set[i].rate * (set[1].hi_start - set[1].lo_start);
    } else if (set[i].style == TRATE) {
      set[i].tilt_stop = set[i].tilt_start * exp(set[i].rate * delt);
    } else if (set[i].style == WIGGLE) {
      double shift = set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].tilt_stop = set[i].tilt_start + shift;

      // compute min/max for WIGGLE = extrema tilt factor will ever reach

      if (set[i].amplitude >= 0.0) {
        if (delt < 0.25 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start;
          set[i].tilt_max = set[i].tilt_start + shift;
        } else if (delt < 0.5 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        } else if (delt < 0.75 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - shift;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        } else {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        }
      } else {
        if (delt < 0.25 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - shift;
          set[i].tilt_max = set[i].tilt_start;
        } else if (delt < 0.5 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start;
        } else if (delt < 0.75 * set[i].tperiod) {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + shift;
        } else {
          set[i].tilt_min = set[i].tilt_start - set[i].amplitude;
          set[i].tilt_max = set[i].tilt_start + set[i].amplitude;
        }
      }
    }
  }

  // if using tilt TRATE, then initial tilt must be non-zero

  for (int i = 3; i < 6; i++)
    if (set[i].style == TRATE && set[i].tilt_start == 0.0)
      error->all(FLERR,"Cannot use fix deform trate on a box with zero tilt");

  // if yz changes and will cause box flip, then xy cannot be changing
  // yz = [3], xy = [5]
  // this is b/c the flips would induce continuous changes in xz
  //   in order to keep the edge vectors of the flipped shape matrix
  //   an integer combination of the edge vectors of the unflipped shape matrix
  // VARIABLE or PRESSURE for yz is error, since no way to calculate if box flip occurs
  // WIGGLE lo/hi flip test is on min/max oscillation limit, not tilt_stop
  // only trigger actual errors if flipflag is set

  if (set[3].style && set[5].style) {
    int flag = 0;
    double lo,hi;
    if (flipflag && (set[3].style == VARIABLE || set[3].style == PRESSURE))
      error->all(FLERR,"Fix deform cannot use yz variable or pressure with xy");
    if (set[3].style == WIGGLE) {
      lo = set[3].tilt_min;
      hi = set[3].tilt_max;
    } else lo = hi = set[3].tilt_stop;
    if (flipflag) {
      if (lo / (set[1].hi_start - set[1].lo_start) < -0.5 ||
          hi / (set[1].hi_start - set[1].lo_start) > 0.5) flag = 1;
      if (set[1].style) {
        if (lo / (set[1].hi_stop - set[1].lo_stop) < -0.5 ||
            hi / (set[1].hi_stop - set[1].lo_stop) > 0.5) flag = 1;
      }
      if (flag)
        error->all(FLERR,"Fix deform is changing yz too much with xy");
    }
  }

  // set domain->h_rate values for use by domain and other fixes/computes
  // initialize all rates to 0.0
  // cannot set here for TRATE,VOLUME,WIGGLE,VARIABLE,PRESSURE since not constant

  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;

  for (int i = 0; i < 3; i++) {
    h_rate[i] = h_ratelo[i] = 0.0;
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == SCALE || set[i].style == VEL ||
        set[i].style == ERATE) {
      double dlo_dt,dhi_dt;
      if (delt != 0.0) {
        dlo_dt = (set[i].lo_stop - set[i].lo_start) / delt;
        dhi_dt = (set[i].hi_stop - set[i].hi_start) / delt;
      } else dlo_dt = dhi_dt = 0.0;
      h_rate[i] = dhi_dt - dlo_dt;
      h_ratelo[i] = dlo_dt;
    }
  }

  for (int i = 3; i < 6; i++) {
    h_rate[i] = 0.0;
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == VEL || set[i].style == ERATE) {
      if (delt != 0.0)
        h_rate[i] = (set[i].tilt_stop - set[i].tilt_start) / delt;
      else h_rate[i] = 0.0;
    }
  }

  // detect if any rigid fixes exist so rigid bodies can be rescaled
  // rfix[] = vector with pointers to each fix rigid

  rfix.clear();

  for (auto ifix : modify->get_fix_list())
    if (ifix->rigid_flag) rfix.push_back(ifix);

  // Find pressure/temp computes if needed

  if (pressure_flag) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Temperature ID for fix deform does not exist");
    temperature = modify->compute[icompute];

    icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all(FLERR,"Pressure ID for fix deform does not exist");
    pressure = modify->compute[icompute];
  }
}

/* ----------------------------------------------------------------------
   compute T,P if needed before integrator starts
------------------------------------------------------------------------- */

void FixDeform::setup(int /*vflag*/)
{
  // trigger virial computation on next timestep
  if (pressure_flag) pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  reset box tilts for flipped config and create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

void FixDeform::pre_exchange()
{
  if (flip == 0) return;

  domain->yz = set[3].tilt_target = set[3].tilt_flip;
  domain->xz = set[4].tilt_target = set[4].tilt_flip;
  domain->xy = set[5].tilt_target = set[5].tilt_flip;
  domain->set_global_box();
  domain->set_local_box();

  domain->image_flip(flipxy,flipxz,flipyz);

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  domain->x2lamda(atom->nlocal);
  irregular->migrate_atoms();
  domain->lamda2x(atom->nlocal);

  flip = 0;
}

/* ---------------------------------------------------------------------- */

void FixDeform::end_of_step()
{
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  // set new box size for strain-based dims

  set_strain();

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

  // Save pressure/strain rate if required

  if (pressure_flag) {
    for (int i = 0; i < 3; i++) {
      set[i].prior_pressure = pressure->vector[i];
      set[i].prior_rate = ((set[i].hi_target - set[i].lo_target) /
                           (domain->boxhi[i] - domain->boxlo[i]) - 1.0)  / update->dt;
    }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);

  // tilt_target can be large positive or large negative value
  // add/subtract box lengths until tilt_target is closest to current value

  if (triclinic) {
    double *h = domain->h;
    for (int i = 3; i < 6; i++) {
      int idenom = 0;
      if (i == 5) idenom = 0;
      else if (i == 4) idenom = 0;
      else if (i == 3) idenom = 1;
      double denom = set[idenom].hi_target - set[idenom].lo_target;

      double current = h[i] / h[idenom];

      while (set[i].tilt_target / denom - current > 0.0)
        set[i].tilt_target -= denom;
      while (set[i].tilt_target / denom - current < 0.0)
        set[i].tilt_target += denom;
      if (fabs(set[i].tilt_target / denom - 1.0 - current) <
          fabs(set[i].tilt_target / denom - current))
        set[i].tilt_target -= denom;
    }
  }

  // if any tilt ratios exceed 0.5, set flip = 1 and compute new tilt values
  // do not flip in x or y if non-periodic (can tilt but not flip)
  //   this is b/c the box length would be changed (dramatically) by flip
  // if yz tilt exceeded, adjust C vector by one B vector
  // if xz tilt exceeded, adjust C vector by one A vector
  // if xy tilt exceeded, adjust B vector by one A vector
  // check yz first since it may change xz, then xz check comes after
  // flip is performed on next timestep, before reneighboring in pre-exchange()

  if (triclinic && flipflag) {
    double xprd = set[0].hi_target - set[0].lo_target;
    double yprd = set[1].hi_target - set[1].lo_target;
    double xprdinv = 1.0 / xprd;
    double yprdinv = 1.0 / yprd;
    if (set[3].tilt_target * yprdinv < -0.5 ||
        set[3].tilt_target * yprdinv > 0.5 ||
        set[4].tilt_target * xprdinv < -0.5 ||
        set[4].tilt_target * xprdinv > 0.5 ||
        set[5].tilt_target * xprdinv < -0.5 ||
        set[5].tilt_target * xprdinv > 0.5) {
      set[3].tilt_flip = set[3].tilt_target;
      set[4].tilt_flip = set[4].tilt_target;
      set[5].tilt_flip = set[5].tilt_target;

      flipxy = flipxz = flipyz = 0;

      if (domain->yperiodic) {
        if (set[3].tilt_flip * yprdinv < -0.5) {
          set[3].tilt_flip += yprd;
          set[4].tilt_flip += set[5].tilt_flip;
          flipyz = 1;
        } else if (set[3].tilt_flip * yprdinv > 0.5) {
          set[3].tilt_flip -= yprd;
          set[4].tilt_flip -= set[5].tilt_flip;
          flipyz = -1;
        }
      }
      if (domain->xperiodic) {
        if (set[4].tilt_flip * xprdinv < -0.5) {
          set[4].tilt_flip += xprd;
          flipxz = 1;
        }
        if (set[4].tilt_flip * xprdinv > 0.5) {
          set[4].tilt_flip -= xprd;
          flipxz = -1;
        }
        if (set[5].tilt_flip * xprdinv < -0.5) {
          set[5].tilt_flip += xprd;
          flipxy = 1;
        }
        if (set[5].tilt_flip * xprdinv > 0.5) {
          set[5].tilt_flip -= xprd;
          flipxy = -1;
        }
      }

      flip = 0;
      if (flipxy || flipxz || flipyz) flip = 1;
      if (flip) next_reneighbor = update->ntimestep + 1;
    }
  }

  // convert atoms and rigid bodies to lamda coords

  if (remapflag == Domain::X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        domain->x2lamda(x[i],x[i]);

    for (auto ifix : rfix)
      ifix->deform(0);
  }

  // reset global and local box to new size/shape
  // only if deform fix is controlling the dimension

  if (set[0].style) {
    domain->boxlo[0] = set[0].lo_target;
    domain->boxhi[0] = set[0].hi_target;
  }
  if (set[1].style) {
    domain->boxlo[1] = set[1].lo_target;
    domain->boxhi[1] = set[1].hi_target;
  }
  if (set[2].style) {
    domain->boxlo[2] = set[2].lo_target;
    domain->boxhi[2] = set[2].hi_target;
  }
  if (triclinic) {
    if (set[3].style) domain->yz = set[3].tilt_target;
    if (set[4].style) domain->xz = set[4].tilt_target;
    if (set[5].style) domain->xy = set[5].tilt_target;
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms and rigid bodies back to box coords

  if (remapflag == Domain::X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        domain->lamda2x(x[i],x[i]);

    for (auto ifix : rfix)
      ifix->deform(1);
  }

  // redo KSpace coeffs since box has changed

  if (kspace_flag) force->kspace->setup();

  // trigger virial computation, if needed, on next timestep

  if (pressure_flag)
    pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   set box size for strain-based dimensions
------------------------------------------------------------------------- */

void FixDeform::set_strain()
{
  // for NONE, target is current box size
  // for TRATE, set target directly based on current time, also set h_rate
  // for WIGGLE, set target directly based on current time, also set h_rate
  // for VARIABLE, set target directly via variable eval, also set h_rate
  // for others except VOLUME, target is linear value between start and stop

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (set[i].style == NONE) {
      set[i].lo_target = domain->boxlo[i];
      set[i].hi_target = domain->boxhi[i];
    } else if (set[i].style == TRATE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      double shift = 0.5 * ((set[i].hi_start - set[i].lo_start) * exp(set[i].rate * delt));
      set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
      set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
      h_rate[i] = set[i].rate * domain->h[i];
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style == WIGGLE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      double shift = 0.5 * set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
      set[i].lo_target = set[i].lo_start - shift;
      set[i].hi_target = set[i].hi_start + shift;
      h_rate[i] = MY_2PI / set[i].tperiod * set[i].amplitude *
        cos(MY_2PI * delt / set[i].tperiod);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style == VARIABLE) {
      double del = input->variable->compute_equal(set[i].hvar);
      set[i].lo_target = set[i].lo_start - 0.5 * del;
      set[i].hi_target = set[i].hi_start + 0.5 * del;
      h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (set[i].style != VOLUME) {
      set[i].lo_target = set[i].lo_start + delta * (set[i].lo_stop - set[i].lo_start);
      set[i].hi_target = set[i].hi_start + delta * (set[i].hi_stop - set[i].hi_start);
    }
  }

  // for triclinic, set new box shape
  // for NONE, target is current tilt
  // for TRATE, set target directly based on current time. also set h_rate
  // for WIGGLE, set target directly based on current time. also set h_rate
  // for VARIABLE, set target directly via variable eval. also set h_rate
  // for other styles, target is linear value between start and stop values

  if (triclinic) {
    double *h = domain->h;

    for (int i = 3; i < 6; i++) {
      if (set[i].style == NONE) {
        if (i == 5) set[i].tilt_target = domain->xy;
        else if (i == 4) set[i].tilt_target = domain->xz;
        else if (i == 3) set[i].tilt_target = domain->yz;
      } else if (set[i].style == TRATE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        set[i].tilt_target = set[i].tilt_start * exp(set[i].rate * delt);
        h_rate[i] = set[i].rate * domain->h[i];
      } else if (set[i].style == WIGGLE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        set[i].tilt_target = set[i].tilt_start +
          set[i].amplitude * sin(MY_2PI * delt / set[i].tperiod);
        h_rate[i] = MY_2PI / set[i].tperiod * set[i].amplitude *
          cos(MY_2PI * delt / set[i].tperiod);
      } else if (set[i].style == VARIABLE) {
        double delta_tilt = input->variable->compute_equal(set[i].hvar);
        set[i].tilt_target = set[i].tilt_start + delta_tilt;
        h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      } else {
        set[i].tilt_target = set[i].tilt_start + delta * (set[i].tilt_stop - set[i].tilt_start);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set box size for pressure-based dimensions
------------------------------------------------------------------------- */

void FixDeform::set_pressure()
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
      if (fabs(set[i].ptarget) > max_h_rate)
        h_rate[i] = max_h_rate * h_rate[i] / fabs(h_rate[i]);

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

void FixDeform::set_volume()
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
        double e3 = (L3 / L3i - 1.0) / dt;
        double p1 = pressure->vector[i];
        double p2 = pressure->vector[fixed];
        double p1i = set[i].prior_pressure;
        double p2i = set[fixed].prior_pressure;

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
            e1 = -e3 / (1 + e3 * dt) * (p2 - p2i) - e2i * (p1 - p2) / (p2 - p2i + (p1 - p1i) / e1i * e2i);
            e2 = (1.0 - (1 + e3 * dt) * (1 + e1 * dt)) / ((1 + e3 * dt) * (1 + e1 * dt) * dt);

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

    set[i].lo_target = 0.5 * (set[i].lo_start + set[i].hi_start) - shift;
    set[i].hi_target = 0.5 * (set[i].lo_start + set[i].hi_start) + shift;
  }
}

/* ----------------------------------------------------------------------
   write Set data to restart file
------------------------------------------------------------------------- */

void FixDeform::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 6 * sizeof(Set);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(set,sizeof(Set),6,fp);
  }
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixDeform::restart(char *buf)
{
  int samestyle = 1;
  Set *set_restart = (Set *) buf;
  for (int i=0; i<6; ++i) {
    // restore data from initial state
    set[i].lo_initial = set_restart[i].lo_initial;
    set[i].hi_initial = set_restart[i].hi_initial;
    set[i].vol_initial = set_restart[i].vol_initial;
    set[i].tilt_initial = set_restart[i].tilt_initial;
    set[i].saved = set_restart[i].saved;
    set[i].prior_rate = set_restart[i].prior_rate;
    set[i].prior_pressure = set_restart[i].prior_pressure;
    // check if style settings are consistent (should do the whole set?)
    if (set[i].style != set_restart[i].style)
      samestyle = 0;
    if (set[i].substyle != set_restart[i].substyle)
      samestyle = 0;
  }
  if (!samestyle)
    error->all(FLERR,"Fix deform settings not consistent with restart");
}

/* ---------------------------------------------------------------------- */

void FixDeform::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix deform command");

  remapflag = Domain::X_REMAP;
  scaleflag = 1;
  flipflag = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform remap", error);
      if (strcmp(arg[iarg+1],"x") == 0) remapflag = Domain::X_REMAP;
      else if (strcmp(arg[iarg+1],"v") == 0) remapflag = Domain::V_REMAP;
      else if (strcmp(arg[iarg+1],"none") == 0) remapflag = Domain::NO_REMAP;
      else error->all(FLERR,"Illegal fix deform remap command: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform units", error);
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix deform units command: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"flip") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform flip", error);
      flipflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform couple", error);
      if (strcmp(arg[iarg+1],"xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1],"xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1],"yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1],"xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1],"none") == 0) pcouple = NOCOUPLE;
      else error->all(FLERR,"Illegal fix fix deform couple command: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max/rate") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform max/rate", error);
      max_h_rate = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (max_h_rate <= 0.0)
        error->all(FLERR,"Maximum strain rate must be a positive, non-zero value");
      iarg += 2;
    } else if (strcmp(arg[iarg],"normalize/pressure") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform normalize/pressure", error);
      normalize_pressure_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vol/balance/p") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix deform vol/balance/p", error);
      vol_balance_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix deform command: {}", arg[iarg]);
  }

  if (dimension == 2)
    if (pcouple == XYZ || pcouple == XZ || pcouple == YZ)
      error->all(FLERR, "Cannot couple Z dimension in fix deform in 2D");
}

/* ----------------------------------------------------------------------
   memory usage of Irregular
------------------------------------------------------------------------- */

double FixDeform::memory_usage()
{
  double bytes = 0.0;
  if (irregular) bytes += irregular->memory_usage();
  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixDeform::modify_param(int narg, char **arg)
{
  if (!pressure_flag) error->all(FLERR,"Cannot modify fix deform without a pressure control");
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify deform", error);
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for deform is not for group all");

    // reset id_temp of pressure to new temperature ID

    icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR,"Pressure ID for fix deform does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    id_press = utils::strdup(arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    return 2;
  }

  return 0;
}
