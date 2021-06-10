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
   Contributing author: Zheng GONG (ENS de Lyon, z.gong@outlook.com)
------------------------------------------------------------------------- */

#include "compute_viscosity_cos.h"

#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeViscosityCos::ComputeViscosityCos(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg) {
  if (narg != 3) error->all(FLERR, "Illegal compute viscosity/cos command");

  scalar_flag = vector_flag = 1;
  size_vector = 7;
  extscalar = 0;
  extvector = -1;
  extlist = new int[7]{1,1,1,1,1,1,0};
  tempflag = 1;
  tempbias = 1;

  maxbias = 0;
  vbiasall = nullptr;

  vector = new double[7];
}

/* ---------------------------------------------------------------------- */

ComputeViscosityCos::~ComputeViscosityCos() {
  if (!copymode) {
    delete[] vector;
    delete[] extlist;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeViscosityCos::setup() {
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeViscosityCos::dof_compute() {
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */
void ComputeViscosityCos::calc_V() {
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone;

  double V_m[2];
  double V_m_local[2] = {0, 0};

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      V_m_local[0] +=
          2 * massone * v[i][0] * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
      V_m_local[1] += massone;
    }

  MPI_Allreduce(V_m_local, V_m, 2, MPI_DOUBLE, MPI_SUM, world);
  V = V_m[0] / V_m[1];
}

/* ---------------------------------------------------------------------- */

double ComputeViscosityCos::compute_scalar() {
  invoked_scalar = update->ntimestep;

  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  double vx_acc;
  double massone;

  calc_V();

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      vx_acc = V * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
      t += ((v[i][0] - vx_acc) * (v[i][0] - vx_acc) + v[i][1] * v[i][1] +
            v[i][2] * v[i][2]) * massone;
    }

  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeViscosityCos::compute_vector() {
  int i;

  invoked_vector = update->ntimestep;

  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vx_acc;

  double massone, t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      vx_acc = V * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
      t[0] += massone * (v[i][0] - vx_acc) * (v[i][0] - vx_acc);
      t[1] += massone * v[i][1] * v[i][1];
      t[2] += massone * v[i][2] * v[i][2];
      t[3] += massone * (v[i][0] - vx_acc) * v[i][1];
      t[4] += massone * (v[i][0] - vx_acc) * v[i][2];
      t[5] += massone * v[i][1] * v[i][2];
    }

  MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
  vector[6] = V;
}


/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeViscosityCos::remove_bias(int i, double *v) {

  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;

  vbias[0] = V * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
  vbias[1] = 0;
  vbias[2] = 0;
  v[0] -= vbias[0];
//  v[1] -= vbias[1];
//  v[2] -= vbias[2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeViscosityCos::remove_bias_thr(int i, double *v, double *b) {
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;

  b[0] = V * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
  b[1] = 0;
  b[2] = 0;
  v[0] -= b[0];
//  v[1] -= b[1];
//  v[2] -= b[2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeViscosityCos::remove_bias_all() {
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vbiasall[i][0] = V * cos(MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
      vbiasall[i][1] = 0;
      vbiasall[i][2] = 0;
      v[i][0] -= vbiasall[i][0];
//      v[i][1] -= vbiasall[i][1];
//      v[i][2] -= vbiasall[i][2];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeViscosityCos::restore_bias(int /* i */, double *v) {
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias_thr()
   assume remove_bias_thr() was previously called with the same buffer b
------------------------------------------------------------------------- */

void ComputeViscosityCos::restore_bias_thr(int /* i */, double *v, double *b) {
  v[0] += b[0];
  v[1] += b[1];
  v[2] += b[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeViscosityCos::restore_bias_all() {
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

