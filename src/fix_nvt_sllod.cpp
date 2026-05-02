// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   LAMMPS development team: developers@lammps.org, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "fix_nvt_sllod.h"

#include "atom.h"
#include "comm.h"
#include "compute_temp_deform.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTSllod::FixNVTSllod(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR, 2, "Temperature control must be used with fix nvt/sllod");
  if (pstat_flag)
    error->all(FLERR, 2, "Pressure control can not be used with fix nvt/sllod");

  // default values

  psllod_flag = 0;
  peculiar_flag = 0;
  integrator = REVERSIBLE;
  bool user_kick = false;
  if (mtchain_default_flag) mtchain = 1;

  // select SLLOD/p-SLLOD/g-SLLOD variant and velocity frame

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"psllod") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} psllod", style), error);
      psllod_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"peculiar") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} peculiar", style), error);
      peculiar_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"kick") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} kick", style), error);
      kick_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      user_kick = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"integrator") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod integrator", error);
      if (strcmp(arg[iarg+1],"reversible") == 0)  integrator = REVERSIBLE;
      else if (strcmp(arg[iarg+1],"legacy") == 0) integrator = LEGACY;
      else error->all(FLERR, "Unknown fix {} integrator argument: {}", style, arg[iarg+1]);
      iarg += 2;
    } else iarg++;
  }

  // default to applying velocity kick in lab frame
  if (!user_kick) kick_flag = !peculiar_flag;

  if (integrator == LEGACY && peculiar_flag == 1)
    error->all(FLERR, "fix {} legacy integrator is not compatible with peculiar=yes", style);

  // create a new compute temp style
  // id = fix-ID + temp

  id_temp = utils::strdup(std::string(id) + "_temp");
  if (peculiar_flag) modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  else modify->add_compute(fmt::format("{} {} temp/deform",id_temp,group->names[igroup]));
  tcomputeflag = 1;
  nondeformbias = 0;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllod::init()
{
  FixNH::init();

  if (!peculiar_flag && !temperature->tempbias)
    error->all(FLERR,"Temperature for fix {} does not have a bias", style);

  nondeformbias = 0;
  if (!utils::strmatch(temperature->style,"^temp/deform")) {
    if (integrator == LEGACY) {
      nondeformbias = 1;
      if (kick_flag)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix {} with peculiar=no and kick=yes requires temperature bias "
                   "to be calculated by compute temp/deform", style);
    } else if (!peculiar_flag) {
      error->all(FLERR, Error::NOLASTLINE, "Fix {} used with lab-frame velocity and non-deform "
                 "temperature bias. For non-deform biases, either set peculiar = yes "
                 "or pass an explicit temp/deform with an extra bias", style);
    }
  }

  // check fix deform remap settings

  auto deform = modify->get_fix_by_style("^deform");
  if (deform.size() < 1)
    error->all(FLERR, Error::NOLASTLINE, "Using fix {} with no fix deform defined", style);

  for (auto &ifix : deform) {
    auto *f = dynamic_cast<FixDeform *>(ifix);
    // not compatible with fix deform. ignore.
    if (!f) continue;
    if ((peculiar_flag && f->remapflag != Domain::NO_REMAP) ||
        (!peculiar_flag && f->remapflag != Domain::V_REMAP))
      error->all(FLERR, Error::NOLASTLINE,
                 "Using fix {} with inconsistent fix {} remap option", style, f->style);

    if (kick_flag) {
      // apply initial kick if velocity stored in lab frame
      // only kick once by default for correct dynamics with multiple run commands
      // make sure fix deform init happens first so h_rate is set
      if (!peculiar_flag) {
        f->init();
        if (comm->me == 0) utils::logmesg(lmp, "fix {} applying velocity profile kick.\n", style);
        auto *f2 = dynamic_cast<ComputeTempDeform *>(temperature);
        if (f2) f2->apply_deform_bias_all();
        kick_flag = 0;
      } else if (comm->me == 0) {
        error->warning(FLERR,"fix {} using peculiar frame velocity. Ignoring kick flag.", style);
      }
    }

    break;
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions with streaming velocity
   also perform sllod update reversibly
-----------------------------------------------------------------------*/

void FixNVTSllod::nve_x()
{
  if (integrator == LEGACY) return FixNH::nve_x();

  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group
  // identical for SLLOD and p-SLLOD
  // velocity treated in peculiar frame relative to sllod streaming for
  //  reversibility, so need to manually account for change in streaming
  //  velocity

  double grad_u[6], xfac[3];
  MathExtra::multiply_shape_shape(domain->h_rate, domain->h_inv, grad_u);
  xfac[0] = exp(grad_u[0]*dthalf);
  xfac[1] = exp(grad_u[1]*dthalf);
  xfac[2] = exp(grad_u[2]*dthalf);
  double vfac[3];
  vfac[0] = exp(-grad_u[0]*dthalf);
  vfac[1] = exp(-grad_u[1]*dthalf);
  vfac[2] = exp(-grad_u[2]*dthalf);

  if (!peculiar_flag)
    dynamic_cast<ComputeTempDeform*>(temperature)->remove_deform_bias_all();

  // fix deform uses the box center as origin for elongation,
  // and the lower corner for shear, so adjust for that
  // to avoid an apparent drift relative to the box and prevent
  // extra atom exchanges between MPI ranks
  double xmid[3];
  for (int i = 0; i < 3; ++i) {
    xmid[i] = (domain->boxhi[i] + domain->boxlo[i])/2.;
  }
  double *xlo = domain->boxlo;

  // propagate boxlo to make second half step reversible
  // xmid does not change
  double ylo2 = xmid[1] + (xlo[1] - xmid[1])*xfac[1]*xfac[1];
  double zlo2 = xmid[2] + (xlo[2] - xmid[2])*xfac[2]*xfac[2];

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      // first half sllod update
      if (psllod_flag) {
        v[i][2] -= dthalf*grad_u[2]*grad_u[2]*x[i][2];
        v[i][1] -= dthalf*grad_u[3]*v[i][2] + dthalf*grad_u[1]*grad_u[1]*x[i][1];
        v[i][0] -= dthalf*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                   + dthalf*grad_u[0]*grad_u[0]*x[i][0];
      } else {
        v[i][1] -= dthalf*grad_u[3]*v[i][2];
        v[i][0] -= dthalf*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
      }
      v[i][0] *= vfac[0];
      v[i][1] *= vfac[1];
      v[i][2] *= vfac[2];

      x[i][1] += dthalf * grad_u[3]*(x[i][2] - xlo[2]);
      x[i][0] += dthalf * (grad_u[5]*(x[i][1] - xlo[1]) + grad_u[4]*(x[i][2] - xlo[2]));
      x[i][0] = xmid[0] + (x[i][0] - xmid[0])*xfac[0];
      x[i][1] = xmid[1] + (x[i][1] - xmid[1])*xfac[1];
      x[i][2] = xmid[2] + (x[i][2] - xmid[2])*xfac[2];

      // nve position update
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // 2nd half sllod update
      x[i][0] = xmid[0] + (x[i][0] - xmid[0])*xfac[0];
      x[i][1] = xmid[1] + (x[i][1] - xmid[1])*xfac[1];
      x[i][2] = xmid[2] + (x[i][2] - xmid[2])*xfac[2];
      x[i][0] += dthalf * (grad_u[5]*(x[i][1] - ylo2) + grad_u[4]*(x[i][2] - zlo2));
      x[i][1] += dthalf * grad_u[3]*(x[i][2] - zlo2);

      // second half sllod velocity step
      // apply here so streaming component matches x when storing in lab frame
      v[i][0] *= vfac[0];
      v[i][1] *= vfac[1];
      v[i][2] *= vfac[2];
      if (psllod_flag) {
        v[i][0] -= dthalf*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2])
                   + dthalf*grad_u[0]*grad_u[0]*x[i][0];
        v[i][1] -= dthalf*grad_u[3]*v[i][2] + dthalf*grad_u[1]*grad_u[1]*x[i][1];
        v[i][2] -= dthalf*grad_u[2]*grad_u[2]*x[i][2];
      } else {
        v[i][0] -= dthalf*(grad_u[5]*v[i][1] + grad_u[4]*v[i][2]);
        v[i][1] -= dthalf*grad_u[3]*v[i][2];
      }
    }
  }

  // x has changed, so can't just call restore_deform_bias_all
  // pass in dtv to account for update to box shape
  if (!peculiar_flag)
    dynamic_cast<ComputeTempDeform*>(temperature)->apply_deform_bias_all(dtv);
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities using legacy method
   NOTE: this style of integration is not time-reversible under mixed
         flows, and neglects the change in streaming velocity caused by
         the position update.
-----------------------------------------------------------------------*/

void FixNVTSllod::nh_v_temp()
{
  if (integrator == REVERSIBLE) return FixNH::nh_v_temp();

  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  if (nondeformbias) temperature->compute_scalar();

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double h_two[6],vdelu[3];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (!psllod_flag) temperature->remove_bias(i,v[i]);
      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      if (psllod_flag) temperature->remove_bias(i,v[i]);
      v[i][0] = v[i][0]*factor_eta - dthalf*vdelu[0];
      v[i][1] = v[i][1]*factor_eta - dthalf*vdelu[1];
      v[i][2] = v[i][2]*factor_eta - dthalf*vdelu[2];
      temperature->restore_bias(i,v[i]);
    }
  }
}

/* ----------------------------------------------------------------------
    calculate the number of data to be packed
------------------------------------------------------------------------- */

int FixNVTSllod::size_restart_global()
{
  return FixNH::size_restart_global() + 1;
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

int FixNVTSllod::pack_restart_data(double *list)
{
  list[0] = kick_flag;
  int n = 1 + FixNH::pack_restart_data(&list[1]);

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixNVTSllod::restart(char *buf)
{
  auto *list = (double *) buf;
  kick_flag = static_cast<int>(list[0]);
  FixNH::restart(buf + sizeof(double));
}

/* ---------------------------------------------------------------------- */

int FixNVTSllod::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"kick") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    kick_flag = utils::logical(FLERR,arg[1],false,lmp);
    return 2;
  }
  return FixNH::modify_param(narg, arg);
}
