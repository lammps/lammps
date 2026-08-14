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
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "compute_temp_deform_eff.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "fix_nh.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { NOBIAS, BIAS };

/* ---------------------------------------------------------------------- */

ComputeTempDeformEff::ComputeTempDeformEff(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), temperature(nullptr), id_temp(nullptr)
{
  tcompute_eff = 0;
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

  if (!atom->electron_flag)
    error->all(FLERR, 2, "Compute {} requires atom style electron", style);

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

void ComputeTempDeformEff::post_constructor() {
  if (tcomputeflag) {
    id_temp = utils::strdup(std::string(id) + "_temp");
    modify->add_compute(fmt::format("{} {} temp/eff", id_temp, group->names[igroup]));
  }
}

/* ---------------------------------------------------------------------- */

ComputeTempDeformEff::~ComputeTempDeformEff()
{
  memory->destroy(vbiasall);

  // delete temperature compute if created by this compute

  if (tcomputeflag) modify->delete_compute(id_temp);
  delete[] id_temp;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::init()
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

  // Flag if internal temperature compute is not an eff compute

  tcompute_eff = 1;
  if (!utils::strmatch(temperature->style, "/eff$"))
    tcompute_eff = 0;

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

void ComputeTempDeformEff::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::dof_compute()
{
  dof = temperature->dof;

  // just include nuclear dof
  // already handled if internal temp compute is an eff style

  if (!tcompute_eff) {
    int *spin = atom->spin;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int one = 0;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && abs(spin[i]) == 1)
        one++;
    int nelectrons;
    MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);

    // Assume 3/2 k T per nucleus

    dof -= domain->dimension * nelectrons;
  }

  // calculate scaling factor to account for different dof used by internal temp compute

  if (dof > 0) tfactor = temperature->dof / dof;
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeformEff::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double *mass = atom->mass;
  if (!mass) {
    scalar = 0;
    return scalar;
  }
  remove_deform_bias_all();
  scalar = temperature->compute_scalar();
  if (dynamic) dof_compute();
  scalar *= tfactor;

  if (!tcompute_eff) {
    int nlocal = atom->nlocal;
    int *spin = atom->spin;
    int *type = atom->type;
    int *mask = atom->mask;
    double *ervel = atom->ervel;
    double mefactor = domain->dimension/4.0;
    double tspin = 0, tspin_all;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && abs(spin[i])==1)
        tspin += mefactor*mass[type[i]]*ervel[i]*ervel[i];

    MPI_Allreduce(&tspin,&tspin_all,1,MPI_DOUBLE,MPI_SUM,world);
    scalar += tspin_all * force->mvv2e / (dof * force->boltz);
  }

  restore_deform_bias_all();
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::compute_vector()
{
  invoked_vector = update->ntimestep;

  remove_deform_bias_all();
  temperature->compute_vector();
  vector = temperature->vector;
  if (dynamic) dof_compute();

  if (!tcompute_eff) {
    double *ervel = atom->ervel;
    double *mass = atom->mass;
    int *spin = atom->spin;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double mefactor = domain->dimension/4.0;
    double massone,tspin[3],tspin_all[3];
    for (int i = 0; i < 3; i++) tspin[i] = 0;

    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && abs(spin[i])==1) {
        massone = mass[type[i]];
        tspin[0] += mefactor * massone * ervel[i]*ervel[i];
        tspin[1] += mefactor * massone * ervel[i]*ervel[i];
        tspin[2] += mefactor * massone * ervel[i]*ervel[i];
      }

    MPI_Allreduce(tspin,tspin_all,3,MPI_DOUBLE,MPI_SUM,world);
    for (int i = 0; i < 3; i++) vector[i] += tspin_all[i] * force->mvv2e;
  }

  restore_deform_bias_all();
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_bias(int i, double *v)
{
  remove_deform_bias(i, v);
  if (which == BIAS) temperature->remove_bias(i, v);
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   NOTE: only removes translational velocity bias from electrons
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_bias_all()
{
  remove_deform_bias_all();
  if (which == BIAS) temperature->remove_bias_all();
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeformEff::restore_bias(int i, double *v)
{
  if (which == BIAS) temperature->restore_bias(i, v);
  restore_deform_bias(i, v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeformEff::restore_bias_all()
{
  if (which == BIAS) temperature->restore_bias_all();
  restore_deform_bias_all();
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I due to deformation only
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_deform_bias(int i, double *v)
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
   remove deform velocity bias from all atoms
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_deform_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall, maxbias, 3, "temp/deform/eff:vbiasall");
  }

  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(atom->x[i], lamda);
      vbiasall[i][0] = h_rate[0] * lamda[0] + h_rate[5] * lamda[1] + h_rate[4] * lamda[2] + h_ratelo[0];
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

void ComputeTempDeformEff::restore_deform_bias(int /*i*/, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in deform velocity bias to all atoms removed by
   remove_deform_bias_all()
   assume remove_deform_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeformEff::restore_deform_bias_all()
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

void ComputeTempDeformEff::apply_deform_bias_all(double dtv)
{
  double ** x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Box may not have been updated yet, so use flow tensor with real coords
  double grad_u[6];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,grad_u);
  double xmid[3];
  xmid[0] = (domain->boxhi[0] + domain->boxlo[0])/2.;
  xmid[1] = (domain->boxhi[1] + domain->boxlo[1])/2.;
  xmid[2] = (domain->boxhi[2] + domain->boxlo[2])/2.;

  // if needed, integrate boxlo to account for box not being updated yet
  // xmid does not change
  double ylo = xmid[1] + (domain->boxlo[1] - xmid[1])*exp(grad_u[1]*dtv);
  double zlo = xmid[2] + (domain->boxlo[2] - xmid[2])*exp(grad_u[2]*dtv);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      apply_deform_bias(v[i], x[i], grad_u, xmid, ylo, zlo);
    }
}

/* ----------------------------------------------------------------------
   add in deform velocity bias to v based on x, grad_u, xmid, ylo and zlo
   does not require remove_deform_bias_all() to be previously called
   box may not have been updated yet, so get flow tensor as input
------------------------------------------------------------------------- */

void ComputeTempDeformEff::apply_deform_bias(double *v, double *x, double *grad_u, double *xmid, double ylo, double zlo)
{
  v[0] += (x[0] - xmid[0]) * grad_u[0] + (x[1] - ylo) * grad_u[5] + (x[2] - zlo) * grad_u[4];
  v[1] += (x[1] - xmid[1]) * grad_u[1] + (x[2] - zlo) * grad_u[3];
  v[2] += (x[2] - xmid[2]) * grad_u[2];
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeformEff::memory_usage()
{
  double bytes = (double)maxbias * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int ComputeTempDeformEff::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR,"compute_modify temp/deform/eff", error);
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
