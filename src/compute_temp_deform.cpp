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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "compute_temp_deform.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_nh.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { NOBIAS, BIAS };

/* ---------------------------------------------------------------------- */

ComputeTempDeform::ComputeTempDeform(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), temperature(nullptr), id_temp(nullptr)
{
  tcomputeflag = 1;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("compute {} temp", style), error);
      delete[] id_temp;
      id_temp = utils::strdup(arg[iarg + 1]);
      tcomputeflag = 0;
      iarg += 2;
    } else error->all(FLERR, iarg, "Unknown compute {} keyword: {}", style, arg[iarg]);
  }

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  maxbias = 0;
  vbiasall = nullptr;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::post_constructor()
{
  if (tcomputeflag) {
    id_temp = utils::strdup(std::string(id) + "_temp");
    modify->add_compute(fmt::format("{} {} temp", id_temp, group->names[igroup]));
  }
}

/* ---------------------------------------------------------------------- */

ComputeTempDeform::~ComputeTempDeform()
{
  if (copymode) return;
  memory->destroy(vbiasall);

  // delete temperature compute if created by this compute

  if (tcomputeflag) modify->delete_compute(id_temp);
  delete[] id_temp;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::init()
{
  // check fix deform remap settings

  auto fixes = modify->get_fix_by_style("^deform");
  if (fixes.size() > 0) {
    auto *f = dynamic_cast<FixDeform *>(fixes[0]);
    if (f && f->remapflag == Domain::X_REMAP && comm->me == 0)
      error->warning(FLERR, "Using compute {} with inconsistent fix deform remap option", style);
  } else {
    if (comm->me == 0)
      error->warning(FLERR, "Using compute {} with no fix deform defined", style);
  }

  // check internal temperature compute

  temperature = modify->get_compute_by_id(id_temp);
  if (!temperature)
    error->all(FLERR, Error::NOLASTLINE,
               "Temperature ID {} for compute {} does not exist", id_temp, style);
  if (temperature->tempflag == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Compute {} temperature ID {} does not compute temperature", style, id_temp);
  if (temperature->igroup != igroup)
    error->all(FLERR, Error::NOLASTLINE,
               "Group of temperature compute with ID {} for compute {} does not match",
               id_temp, style);

  // avoid possibility of self-referential loop

  if (utils::strmatch(temperature->style, "^temp/deform"))
    error->all(FLERR, Error::NOLASTLINE,
               "Compute {} internal temperature compute cannot be of style temp/deform", style);

  if (temperature->tempbias)
    which = BIAS;
  else
    which = NOBIAS;

  // make sure internal temperature compute is called first

  temperature->init();
  temperature->setup();

  vector = temperature->vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::dof_compute()
{
  dof = temperature->dof;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeform::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  remove_deform_bias_all();
  scalar = temperature->compute_scalar();
  if (dynamic) dof_compute();
  restore_deform_bias_all();

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeform::compute_vector()
{
  invoked_vector = update->ntimestep;

  remove_deform_bias_all();
  temperature->compute_vector();
  vector = temperature->vector;
  if (dynamic) dof_compute();
  restore_deform_bias_all();
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_bias(int i, double *v)
{
  remove_deform_bias(i, v);
  if (which == BIAS) temperature->remove_bias(i, v);
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_bias_thr(int i, double *v, double *b)
{
  remove_deform_bias_thr(i, v, b);
  if (which == BIAS) temperature->remove_bias_thr(i, v, b);
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_bias_all()
{
  remove_deform_bias_all();
  if (which == BIAS) temperature->remove_bias_all();
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_bias(int i, double *v)
{
  if (which == BIAS) temperature->restore_bias(i, v);
  restore_deform_bias(i, v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_bias_thr(int i, double *v, double *b)
{
  if (which == BIAS) temperature->restore_bias_thr(i, v, b);
  restore_deform_bias_thr(i, v, b);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_bias_all()
{
  if (which == BIAS) temperature->restore_bias_all();
  restore_deform_bias_all();
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I due to deformation only
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_deform_bias(int i, double *v)
{
  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  domain->x2lamda(atom->x[i], lamda);
  vbias[0] = h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
  vbias[1] = h_rate[1] * lamda[1] + h_rate[3] * lamda[2] + h_ratelo[1];
  vbias[2] = h_rate[2] * lamda[2] + h_ratelo[2];
  v[0] -= vbias[0];
  v[1] -= vbias[1];
  v[2] -= vbias[2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I due to deformation only
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_deform_bias_thr(int i, double *v, double *b)
{
  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  domain->x2lamda(atom->x[i], lamda);
  if (which == NOBIAS) {
    b[0] = h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
    b[1] = h_rate[1] * lamda[1] + h_rate[3] * lamda[2] + h_ratelo[1];
    b[2] = h_rate[2] * lamda[2] + h_ratelo[2];
    v[0] -= b[0];
    v[1] -= b[1];
    v[2] -= b[2];
  } else {
    // b needed by internal temperature compute, so just re-calculate deform bias when restoring
    v[0] -= h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
    v[1] -= h_rate[1] * lamda[1] + h_rate[3] * lamda[2] + h_ratelo[1];
    v[2] -= h_rate[2] * lamda[2] + h_ratelo[2];
  }
}

/* ----------------------------------------------------------------------
   remove deform velocity bias from all atoms
------------------------------------------------------------------------- */

void ComputeTempDeform::remove_deform_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall, maxbias, 3, "temp/deform:vbiasall");
  }

  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(atom->x[i], lamda);
      vbiasall[i][0] =
          h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
      vbiasall[i][1] = h_rate[1] * lamda[1] + h_rate[3] * lamda[2] + h_ratelo[1];
      vbiasall[i][2] = h_rate[2] * lamda[2] + h_ratelo[2];
      v[i][0] -= vbiasall[i][0];
      v[i][1] -= vbiasall[i][1];
      v[i][2] -= vbiasall[i][2];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_deform_bias()
   assume remove_deform_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_deform_bias(int /*i*/, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in deform velocity bias to atom I removed by
   remove_deform_bias_thr()
   assume remove_deform_bias_thr() was previously called with the same
   buffer b
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_deform_bias_thr(int i, double *v, double *b)
{
  if (which == NOBIAS) {
    v[0] += b[0];
    v[1] += b[1];
    v[2] += b[2];
  } else {
    double lamda[3];
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    domain->x2lamda(atom->x[i], lamda);
    v[0] += h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
    v[1] += h_rate[1] * lamda[1] + h_rate[3] * lamda[2] + h_ratelo[1];
    v[2] += h_rate[2] * lamda[2] + h_ratelo[2];
  }
}

/* ----------------------------------------------------------------------
   add back in deform velocity bias to all atoms removed by
   remove_deform_bias_all()
   assume remove_deform_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeform::restore_deform_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] += vbiasall[i][0];
      v[i][1] += vbiasall[i][1];
      v[i][2] += vbiasall[i][2];
    }
}

/* ----------------------------------------------------------------------
   add in deform velocity bias to all atoms based on x
   does not require remove_deform_bias_all() to be previously called
   approximately propagate boxlo by dtv for velocity calculation since
   shear velocity is relative to lower corner
------------------------------------------------------------------------- */

void ComputeTempDeform::apply_deform_bias_all(double dtv)
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Box may not have been updated yet, so use flow tensor with real coords
  double grad_u[6];
  MathExtra::multiply_shape_shape(domain->h_rate, domain->h_inv, grad_u);
  double xmid[3];
  xmid[0] = (domain->boxhi[0] + domain->boxlo[0]) / 2.;
  xmid[1] = (domain->boxhi[1] + domain->boxlo[1]) / 2.;
  xmid[2] = (domain->boxhi[2] + domain->boxlo[2]) / 2.;

  // if needed, integrate boxlo to account for box not being updated yet
  // xmid does not change
  double ylo = xmid[1] + (domain->boxlo[1] - xmid[1]) * exp(grad_u[1] * dtv);
  double zlo = xmid[2] + (domain->boxlo[2] - xmid[2]) * exp(grad_u[2] * dtv);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) { apply_deform_bias(v[i], x[i], grad_u, xmid, ylo, zlo); }
}

/* ----------------------------------------------------------------------
   add in deform velocity bias to v based on x, grad_u, xmid, ylo and zlo
   does not require remove_deform_bias_all() to be previously called
   box may not have been updated yet, so get flow tensor as input
------------------------------------------------------------------------- */

void ComputeTempDeform::apply_deform_bias(double *v, double *x, double *grad_u, double *xmid,
                                          double ylo, double zlo)
{
  v[0] += (x[0] - xmid[0]) * grad_u[0] + (x[1] - ylo) * grad_u[5] + (x[2] - zlo) * grad_u[4];
  v[1] += (x[1] - xmid[1]) * grad_u[1] + (x[2] - zlo) * grad_u[3];
  v[2] += (x[2] - xmid[2]) * grad_u[2];
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeform::memory_usage()
{
  double bytes = 3 * maxbias * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int ComputeTempDeform::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify temp/deform", error);
    if (tcomputeflag) modify->delete_compute(id_temp);
    delete[] id_temp;
    tcomputeflag = 0;
    id_temp = utils::strdup(arg[1]);
    return 2;
  } else if (strcmp(arg[0], "extra/dof") == 0) {
    // Can't set extra/dof of internal temp compute directly,
    // so pass through the modify call
    temperature->modify_params(MIN(narg, 2), arg);
  } else if (strcmp(arg[0], "dynamic/dof") == 0) {
    // Can't set dynamic_user flag of internal temp compute directly,
    // so pass through the modify call
    temperature->modify_params(MIN(narg, 2), arg);
  }
  return 0;
}
