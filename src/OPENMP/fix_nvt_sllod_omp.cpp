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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_nvt_sllod_omp.h"

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

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using dbl3_t = struct { double x,y,z; };

/* ---------------------------------------------------------------------- */

FixNVTSllodOMP::FixNVTSllodOMP(LAMMPS *lmp, int narg, char **arg) :
  FixNHOMP(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR, 2, "Temperature control must be used with fix nvt/sllod/omp");
  if (pstat_flag)
    error->all(FLERR, 2, "Pressure control can not be used with fix nvt/sllod/omp");

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
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod/omp psllod", error);
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

void FixNVTSllodOMP::init()
{
  FixNHOMP::init();

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

void FixNVTSllodOMP::nve_x()
{
  if (integrator == LEGACY) return FixNHOMP::nve_x();

  auto * _noalias const v = (dbl3_t *) atom->v[0];
  auto * _noalias const x = (dbl3_t *) atom->x[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  // x update by full step only for atoms in group
  // identical for SLLOD and p-SLLOD
  // velocity treated in peculiar frame relative to sllod streaming for
  //  reversibility, so need to manually account for change in streaming
  //  velocity

  double grad_u[6];
  MathExtra::multiply_shape_shape(domain->h_rate, domain->h_inv, grad_u);
  const dbl3_t xfac = {
    exp(grad_u[0]*dthalf),
    exp(grad_u[1]*dthalf),
    exp(grad_u[2]*dthalf)
  };
  const dbl3_t vfac = {
    exp(-grad_u[0]*dthalf),
    exp(-grad_u[1]*dthalf),
    exp(-grad_u[2]*dthalf)
  };

  // fix deform uses the box center as origin for elongation,
  // and the lower corner for shear, so adjust for that
  // to avoid an apparent drift relative to the box and prevent
  // extra atom exchanges between MPI ranks
  dbl3_t xmid;
  xmid.x = (domain->boxhi[0] + domain->boxlo[0])/2.;
  xmid.y = (domain->boxhi[1] + domain->boxlo[1])/2.;
  xmid.z = (domain->boxhi[2] + domain->boxlo[2])/2.;
  const double * _noalias const xlo = domain->boxlo;

  // propagate boxlo to make second half step reversible
  // xmid does not change
  const double ylo2 = xmid.y + (xlo[1] - xmid.y)*xfac.y*xfac.y;
  const double zlo2 = xmid.z + (xlo[2] - xmid.z)*xfac.z*xfac.z;

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE LMP_SHARED(grad_u) schedule(static)
#endif
  for (int i = 0; i < nlocal; ++i) {
    double buf[3];
    if (mask[i] & groupbit) {
      if (!peculiar_flag)
        dynamic_cast<ComputeTempDeform*>(temperature)->remove_deform_bias_thr(i,&v[i].x,buf);

      // first half sllod update
      if (psllod_flag) {
        v[i].z -= dthalf*grad_u[2]*grad_u[2]*x[i].z;
        v[i].y -= dthalf*grad_u[3]*v[i].z + dthalf*grad_u[1]*grad_u[1]*x[i].y;
        v[i].x -= dthalf*(grad_u[5]*v[i].y + grad_u[4]*v[i].z)
                   + dthalf*grad_u[0]*grad_u[0]*x[i].x;
      } else {
        v[i].y -= dthalf*grad_u[3]*v[i].z;
        v[i].x -= dthalf*(grad_u[5]*v[i].y + grad_u[4]*v[i].z);
      }
      v[i].x *= vfac.x;
      v[i].y *= vfac.y;
      v[i].z *= vfac.z;

      x[i].y += dthalf * grad_u[3]*(x[i].z - xlo[2]);
      x[i].x += dthalf * (grad_u[5]*(x[i].y - xlo[1]) + grad_u[4]*(x[i].z - xlo[2]));
      x[i].x = xmid.x + (x[i].x - xmid.x)*xfac.x;
      x[i].y = xmid.y + (x[i].y - xmid.y)*xfac.y;
      x[i].z = xmid.z + (x[i].z - xmid.z)*xfac.z;

      // nve position update
      x[i].x += dtv * v[i].x;
      x[i].y += dtv * v[i].y;
      x[i].z += dtv * v[i].z;

      // 2nd half sllod update
      x[i].x = xmid.x + (x[i].x - xmid.x)*xfac.x;
      x[i].y = xmid.y + (x[i].y - xmid.y)*xfac.y;
      x[i].z = xmid.z + (x[i].z - xmid.z)*xfac.z;
      x[i].x += dthalf * (grad_u[5]*(x[i].y - ylo2) + grad_u[4]*(x[i].z - zlo2));
      x[i].y += dthalf * grad_u[3]*(x[i].z - zlo2);

      // second half sllod velocity step
      // apply here so streaming component matches x when storing in lab frame
      v[i].x *= vfac.x;
      v[i].y *= vfac.y;
      v[i].z *= vfac.z;
      if (psllod_flag) {
        v[i].x -= dthalf*(grad_u[5]*v[i].y + grad_u[4]*v[i].z)
                   + dthalf*grad_u[0]*grad_u[0]*x[i].x;
        v[i].y -= dthalf*grad_u[3]*v[i].z + dthalf*grad_u[1]*grad_u[1]*x[i].y;
        v[i].z -= dthalf*grad_u[2]*grad_u[2]*x[i].z;
      } else {
        v[i].x -= dthalf*(grad_u[5]*v[i].y + grad_u[4]*v[i].z);
        v[i].y -= dthalf*grad_u[3]*v[i].z;
      }

      if (!peculiar_flag)
        dynamic_cast<ComputeTempDeform*>(temperature)->apply_deform_bias(&v[i].x, &x[i].x, grad_u, &xmid.x, ylo2, zlo2);
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities using legacy method
   NOTE: this style of integration is not time-reversible under mixed
         flows, and neglects the change in streaming velocity caused by
         the position update.
-----------------------------------------------------------------------*/


void FixNVTSllodOMP::nh_v_temp()
{
  if (integrator == REVERSIBLE) return FixNHOMP::nh_v_temp();

  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  auto * _noalias const v = (dbl3_t *) atom->v[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  if (nondeformbias) temperature->compute_scalar();

  double h_two[6];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE LMP_SHARED(h_two) schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    double vdelu0,vdelu1,vdelu2,buf[3];
    if (mask[i] & groupbit) {
      if (!psllod_flag) temperature->remove_bias_thr(i,&v[i].x,buf);
      vdelu0 = h_two[0]*v[i].x + h_two[5]*v[i].y + h_two[4]*v[i].z;
      vdelu1 = h_two[1]*v[i].y + h_two[3]*v[i].z;
      vdelu2 = h_two[2]*v[i].z;
      if (psllod_flag) temperature->remove_bias_thr(i,&v[i].x,buf);
      v[i].x = v[i].x*factor_eta - dthalf*vdelu0;
      v[i].y = v[i].y*factor_eta - dthalf*vdelu1;
      v[i].z = v[i].z*factor_eta - dthalf*vdelu2;
      temperature->restore_bias_thr(i,&v[i].x,buf);
    }
  }
}

/* ----------------------------------------------------------------------
    calculate the number of data to be packed
------------------------------------------------------------------------- */

int FixNVTSllodOMP::size_restart_global()
{
  return FixNHOMP::size_restart_global() + 1;
}

/* ----------------------------------------------------------------------
   pack restart data
------------------------------------------------------------------------- */

int FixNVTSllodOMP::pack_restart_data(double *list)
{
  list[0] = kick_flag;
  int n = 1 + FixNHOMP::pack_restart_data(&list[1]);

  return n;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixNVTSllodOMP::restart(char *buf)
{
  auto *list = (double *) buf;
  kick_flag = static_cast<int>(list[0]);
  FixNHOMP::restart(buf + sizeof(double));
}

/* ---------------------------------------------------------------------- */

int FixNVTSllodOMP::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"kick") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    kick_flag = utils::logical(FLERR,arg[1],false,lmp);
    return 2;
  }
  return FixNHOMP::modify_param(narg, arg);
}
