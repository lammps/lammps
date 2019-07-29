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

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   based on ComputeTempAsphere
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstring>
#include "compute_temp_body.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ROTATE,ALL};

/* ---------------------------------------------------------------------- */

ComputeTempBody::ComputeTempBody(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_bias(NULL), tbias(NULL), avec(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp/body command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  tempbias = 0;
  id_bias = NULL;
  mode = ALL;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/body command");
      tempbias = 1;
      int n = strlen(arg[iarg+1]) + 1;
      id_bias = new char[n];
      strcpy(id_bias,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/body command");
      if (strcmp(arg[iarg+1],"rotate") == 0) mode = ROTATE;
      else if (strcmp(arg[iarg+1],"all") == 0) mode = ALL;
      else error->all(FLERR,"Illegal compute temp/body command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/body command");
  }

  vector = new double[size_vector];

}

/* ---------------------------------------------------------------------- */

ComputeTempBody::~ComputeTempBody()
{
  delete [] id_bias;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempBody::init()
{
  // error check

  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec)
    error->all(FLERR,"Compute temp/body requires atom style body");

  // check that all particles are finite-size, no point particles allowed

  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (body[i] < 0)
        error->one(FLERR,"Compute temp/body requires bodies");

  if (tempbias) {
    int i = modify->find_compute(id_bias);
    if (i < 0)
      error->all(FLERR,"Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all(FLERR,"Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all(FLERR,"Bias compute does not calculate a velocity bias");
    if (tbias->igroup != igroup)
      error->all(FLERR,"Bias compute group does not match compute group");
    if (strcmp(tbias->style,"temp/region") == 0) tempbias = 2;
    else tempbias = 1;

    // init and setup bias compute because
    // this compute's setup()->dof_compute() may be called first

    tbias->init();
    tbias->setup();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTempBody::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempBody::dof_compute()
{
  adjust_dof_fix();

  // 6 dof for 3d, 3 dof for 2d
  // which dof are included also depends on mode
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  natoms_temp = group->count(igroup);
  int nper;
  if (domain->dimension == 3) {
    if (mode == ALL) nper = 6;
    else nper = 3;
  } else {
    if (mode == ALL) nper = 3;
    else nper = 1;
  }
  dof = nper*natoms_temp;

  // additional adjustments to dof

  if (tempbias == 1) {
    if (mode == ALL) dof -= tbias->dof_remove(-1) * natoms_temp;

  } else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    tbias->dof_remove_pre();

    int count = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (tbias->dof_remove(i)) count++;
    int count_all;
    MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
    dof -= nper*count_all;
  }

  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempBody::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  AtomVecBody::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *inertia,*quat;
  double wbody[3];
  double rot[3][3];

  // sum translational and rotational energy for each particle

  double t = 0.0;

  if (mode == ALL) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];

        // principal moments of inertia

        inertia = bonus[body[i]].inertia;
        quat = bonus[body[i]].quat;

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        if (inertia[0] == 0.0) wbody[0] = 0.0;
        else wbody[0] /= inertia[0];
        if (inertia[1] == 0.0) wbody[1] = 0.0;
        else wbody[1] /= inertia[1];
        if (inertia[2] == 0.0) wbody[2] = 0.0;
        else wbody[2] /= inertia[2];

        t += inertia[0]*wbody[0]*wbody[0] +
          inertia[1]*wbody[1]*wbody[1] + inertia[2]*wbody[2]*wbody[2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        // principal moments of inertia

        inertia = bonus[body[i]].inertia;
        quat = bonus[body[i]].quat;

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        if (inertia[0] == 0.0) wbody[0] = 0.0;
        else wbody[0] /= inertia[0];
        if (inertia[1] == 0.0) wbody[1] = 0.0;
        else wbody[1] /= inertia[1];
        if (inertia[2] == 0.0) wbody[2] = 0.0;
        else wbody[2] /= inertia[2];

        t += inertia[0]*wbody[0]*wbody[0] +
          inertia[1]*wbody[1]*wbody[1] + inertia[2]*wbody[2]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempBody::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  AtomVecBody::Bonus *bonus = avec->bonus;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *inertia,*quat;
  double wbody[3],t[6];
  double rot[3][3];
  double massone;

  // sum translational and rotational energy for each particle

  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (mode == ALL) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        t[0] += massone * v[i][0]*v[i][0];
        t[1] += massone * v[i][1]*v[i][1];
        t[2] += massone * v[i][2]*v[i][2];
        t[3] += massone * v[i][0]*v[i][1];
        t[4] += massone * v[i][0]*v[i][2];
        t[5] += massone * v[i][1]*v[i][2];

        // principal moments of inertia

        inertia = bonus[body[i]].inertia;
        quat = bonus[body[i]].quat;

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        if (inertia[0] == 0.0) wbody[0] = 0.0;
        else wbody[0] /= inertia[0];
        if (inertia[1] == 0.0) wbody[1] = 0.0;
        else wbody[1] /= inertia[1];
        if (inertia[2] == 0.0) wbody[2] = 0.0;
        else wbody[2] /= inertia[2];

        // rotational kinetic energy

        t[0] += inertia[0]*wbody[0]*wbody[0];
        t[1] += inertia[1]*wbody[1]*wbody[1];
        t[2] += inertia[2]*wbody[2]*wbody[2];
        t[3] += inertia[0]*wbody[0]*wbody[1];
        t[4] += inertia[1]*wbody[0]*wbody[2];
        t[5] += inertia[2]*wbody[1]*wbody[2];
      }

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        // principal moments of inertia

        inertia = bonus[body[i]].inertia;
        quat = bonus[body[i]].quat;
        massone = rmass[i];

        // wbody = angular velocity in body frame

        MathExtra::quat_to_mat(quat,rot);
        MathExtra::transpose_matvec(rot,angmom[i],wbody);
        if (inertia[0] == 0.0) wbody[0] = 0.0;
        else wbody[0] /= inertia[0];
        if (inertia[1] == 0.0) wbody[1] = 0.0;
        else wbody[1] /= inertia[1];
        if (inertia[2] == 0.0) wbody[2] = 0.0;
        else wbody[2] /= inertia[2];

        // rotational kinetic energy

        t[0] += inertia[0]*wbody[0]*wbody[0];
        t[1] += inertia[1]*wbody[1]*wbody[1];
        t[2] += inertia[2]*wbody[2]*wbody[2];
        t[3] += inertia[0]*wbody[0]*wbody[1];
        t[4] += inertia[1]*wbody[0]*wbody[2];
        t[5] += inertia[2]*wbody[1]*wbody[2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempBody::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempBody::restore_bias(int i, double *v)
{
  if (tbias) tbias->restore_bias(i,v);
}
