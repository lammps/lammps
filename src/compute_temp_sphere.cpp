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

#include "compute_temp_sphere.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ROTATE,ALL};

#define INERTIA 0.4          // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeTempSphere::ComputeTempSphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_bias(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute temp/sphere command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  tempbias = 0;
  mode = ALL;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"bias") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/sphere command");
      tempbias = 1;
      int n = strlen(arg[iarg+1]) + 1;
      id_bias = new char[n];
      strcpy(id_bias,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dof") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/sphere command");
      if (strcmp(arg[iarg+1],"rotate") == 0) mode = ROTATE;
      else if (strcmp(arg[iarg+1],"all") == 0) mode = ALL;
      else error->all(FLERR,"Illegal compute temp/sphere command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/sphere command");
  }

  // when computing only the rotational temperature,
  // do not remove DOFs for translation as set by default

  if (mode == ROTATE) extra_dof = 0;

  vector = new double[size_vector];

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Compute temp/sphere requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

ComputeTempSphere::~ComputeTempSphere()
{
  delete [] id_bias;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::init()
{
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

void ComputeTempSphere::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::dof_compute()
{
  int count,count_all;

  adjust_dof_fix();
  natoms_temp = group->count(igroup);

  // 6 or 3 dof for extended/point particles for 3d
  // 3 or 2 dof for extended/point particles for 2d
  // which dof are included also depends on mode
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  count = 0;
  if (domain->dimension == 3) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (radius[i] == 0.0) {
          if (mode == ALL) count += 3;
        } else {
          if (mode == ALL) count += 6;
          else count += 3;
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (radius[i] == 0.0) {
          if (mode == ALL) count += 2;
        } else {
          if (mode == ALL) count += 3;
          else count += 1;
        }
      }
  }

  MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
  dof = count_all;

  // additional adjustments to dof

  if (tempbias == 1) {
    if (mode == ALL) dof -= tbias->dof_remove(-1) * natoms_temp;

  } else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    tbias->dof_remove_pre();

    count = 0;
    if (domain->dimension == 3) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (tbias->dof_remove(i)) {
            if (radius[i] == 0.0) {
              if (mode == ALL) count += 3;
            } else {
              if (mode == ALL) count += 6;
              else count += 3;
            }
          }
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (tbias->dof_remove(i)) {
            if (radius[i] == 0.0) {
              if (mode == ALL) count += 2;
            } else {
              if (mode == ALL) count += 3;
              else count += 1;
            }
          }
        }
    }

    MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
    dof -= count_all;
  }

  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempSphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // point particles will not contribute rotation due to radius = 0

  double t = 0.0;

  if (mode == ALL) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
        t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
              omega[i][2]*omega[i][2]) * INERTIA*rmass[i]*radius[i]*radius[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
              omega[i][2]*omega[i][2]) * INERTIA*rmass[i]*radius[i]*radius[i];
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

void ComputeTempSphere::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // point particles will not contribute rotation due to radius = 0

  double massone,inertiaone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  if (mode == ALL) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        t[0] += massone * v[i][0]*v[i][0];
        t[1] += massone * v[i][1]*v[i][1];
        t[2] += massone * v[i][2]*v[i][2];
        t[3] += massone * v[i][0]*v[i][1];
        t[4] += massone * v[i][0]*v[i][2];
        t[5] += massone * v[i][1]*v[i][2];

        inertiaone = INERTIA*rmass[i]*radius[i]*radius[i];
        t[0] += inertiaone * omega[i][0]*omega[i][0];
        t[1] += inertiaone * omega[i][1]*omega[i][1];
        t[2] += inertiaone * omega[i][2]*omega[i][2];
        t[3] += inertiaone * omega[i][0]*omega[i][1];
        t[4] += inertiaone * omega[i][0]*omega[i][2];
        t[5] += inertiaone * omega[i][1]*omega[i][2];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        inertiaone = INERTIA*rmass[i]*radius[i]*radius[i];
        t[0] += inertiaone * omega[i][0]*omega[i][0];
        t[1] += inertiaone * omega[i][1]*omega[i][1];
        t[2] += inertiaone * omega[i][2]*omega[i][2];
        t[3] += inertiaone * omega[i][0]*omega[i][1];
        t[4] += inertiaone * omega[i][0]*omega[i][2];
        t[5] += inertiaone * omega[i][1]*omega[i][2];
      }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias(int i, double *v)
{
  tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias_thr(int i, double *v, double *b)
{
  tbias->remove_bias_thr(i,v,b);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias(int i, double *v)
{
  tbias->restore_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias_thr(int i, double *v, double *b)
{
  tbias->restore_bias_thr(i,v,b);
}
