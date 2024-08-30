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

#include "fix_press_berendsen.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NOBIAS, BIAS };
enum { NONE, XYZ, XY, YZ, XZ };
enum { ISO, ANISO };

/* ---------------------------------------------------------------------- */

FixPressBerendsen::FixPressBerendsen(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_temp(nullptr), id_press(nullptr), tflag(0), pflag(0)
{
  if (narg < 5) error->all(FLERR, "Illegal fix press/berendsen command");

  // Berendsen barostat applied every step

  nevery = 1;

  // default values

  pcouple = NONE;
  bulkmodulus = 10.0;
  allremap = 1;

  for (int i = 0; i < 3; i++) {
    p_start[i] = p_stop[i] = p_period[i] = 0.0;
    p_flag[i] = 0;
    p_period[i] = 0.0;
  }

  // process keywords

  dimension = domain->dimension;

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "iso") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      pcouple = XYZ;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg], "aniso") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      pcouple = NONE;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg], "x") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "y") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[1] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[1] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "z") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[2] = 1;
      iarg += 4;
      if (dimension == 2) error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");

    } else if (strcmp(arg[iarg], "couple") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      if (strcmp(arg[iarg + 1], "xyz") == 0)
        pcouple = XYZ;
      else if (strcmp(arg[iarg + 1], "xy") == 0)
        pcouple = XY;
      else if (strcmp(arg[iarg + 1], "yz") == 0)
        pcouple = YZ;
      else if (strcmp(arg[iarg + 1], "xz") == 0)
        pcouple = XZ;
      else if (strcmp(arg[iarg + 1], "none") == 0)
        pcouple = NONE;
      else
        error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;

    } else if (strcmp(arg[iarg], "modulus") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      bulkmodulus = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (bulkmodulus <= 0.0) error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dilate") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix press/berendsen command");
      if (strcmp(arg[iarg + 1], "all") == 0)
        allremap = 1;
      else if (strcmp(arg[iarg + 1], "partial") == 0)
        allremap = 0;
      else
        error->all(FLERR, "Illegal fix press/berendsen command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix press/berendsen command");
  }

  if (allremap == 0) restart_pbc = 1;

  // error checks

  if (dimension == 2 && p_flag[2])
    error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR, "Invalid fix press/berendsen for a 2d simulation");

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR, "Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR, "Cannot use fix press/berendsen on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR, "Cannot use fix press/berendsen on a non-periodic dimension");

  if (pcouple == XYZ && dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] || p_stop[0] != p_stop[1] ||
       p_stop[0] != p_stop[2] || p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XYZ && dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] || p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] || p_period[0] != p_period[1]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] || p_period[1] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] || p_period[0] != p_period[2]))
    error->all(FLERR, "Invalid fix press/berendsen pressure settings");

  if ((p_flag[0] && p_period[0] <= 0.0) || (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0))
    error->all(FLERR, "Fix press/berendsen damping parameters must be > 0.0");

  if (p_flag[0]) box_change |= BOX_CHANGE_X;
  if (p_flag[1]) box_change |= BOX_CHANGE_Y;
  if (p_flag[2]) box_change |= BOX_CHANGE_Z;

  // pstyle = ISO if XYZ coupling or XY coupling in 2d -> 1 dof
  // else pstyle = ANISO -> 3 dof

  if (pcouple == XYZ || (dimension == 2 && pcouple == XY))
    pstyle = ISO;
  else
    pstyle = ANISO;

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp", id_temp));
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure {}", id_press, id_temp));
  pflag = 1;
}

/* ---------------------------------------------------------------------- */

FixPressBerendsen::~FixPressBerendsen()
{
  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete[] id_temp;
  delete[] id_press;
}

/* ---------------------------------------------------------------------- */

int FixPressBerendsen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::init()
{
  if (domain->triclinic) error->all(FLERR, "Cannot use fix press/berendsen with triclinic box");

  // ensure no conflict with fix deform

  for (const auto &ifix : modify->get_fix_by_style("^deform")) {
    int *dimflag = static_cast<FixDeform *>(ifix)->dimflag;
    if (!dimflag) continue;
    if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) || (p_flag[2] && dimflag[2]))
      error->all(FLERR,
                 "Cannot use fix press/berendsen and "
                 "fix deform on same component of stress tensor");
  }

  // set temperature and pressure ptrs

  temperature = modify->get_compute_by_id(id_temp);
  if (!temperature)
    error->all(FLERR, "Temperature compute ID {} for fix press/berendsen does not exist", id_temp);

  if (temperature->tempbias)
    which = BIAS;
  else
    which = NOBIAS;

  pressure = modify->get_compute_by_id(id_press);
  if (!pressure)
    error->all(FLERR, "Pressure compute ID {} for fix press/berendsen does not exist", id_press);

  // Kspace setting

  if (force->kspace)
    kspace_flag = 1;
  else
    kspace_flag = 0;

  // detect if any rigid fixes exist so rigid bodies move when box is remapped

  rfix.clear();
  for (auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag) rfix.push_back(ifix);
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */

void FixPressBerendsen::setup(int /*vflag*/)
{
  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::end_of_step()
{
  // compute new T,P

  if (pstyle == ISO) {
    temperature->compute_scalar();
    pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);
      dilation[i] = pow(1.0 - update->dt / p_period[i] * (p_target[i] - p_current[i]) / bulkmodulus,
                        1.0 / 3.0);
    }
  }

  // remap simulation box and atoms
  // redo KSpace coeffs since volume has changed

  remap();
  if (kspace_flag) force->kspace->setup();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPressBerendsen::couple()
{
  double *tensor = pressure->vector;

  if (pstyle == ISO)
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  else if (pcouple == XYZ) {
    double ave = 1.0 / 3.0 * (tensor[0] + tensor[1] + tensor[2]);
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
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixPressBerendsen::remap()
{
  int i;
  double oldlo, oldhi, ctr;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap)
    domain->x2lamda(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->x2lamda(x[i], x[i]);
  }

  for (auto &ifix : rfix) ifix->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo - ctr) * dilation[i] + ctr;
      domain->boxhi[i] = (oldhi - ctr) * dilation[i] + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap)
    domain->lamda2x(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->lamda2x(x[i], x[i]);
  }

  for (auto &ifix : rfix) ifix->deform(1);
}

/* ---------------------------------------------------------------------- */

int FixPressBerendsen::modify_param(int narg, char **arg)
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
      error->warning(FLERR, "Temperature compute {} for fix {} is not for group all: {}", arg[1],
                     style, group->names[temperature->igroup]);

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
