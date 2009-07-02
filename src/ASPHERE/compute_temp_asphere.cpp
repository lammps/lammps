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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "compute_temp_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::ComputeTempAsphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3 && narg != 4)
    error->all("Illegal compute temp/asphere command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  tempbias = 0;
  id_bias = NULL;
  if (narg == 4) {
    tempbias = 1;
    int n = strlen(arg[3]) + 1;
    id_bias = new char[n];
    strcpy(id_bias,arg[3]);
  }

  vector = new double[6];
  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_temp_sphere:inertia");

  // error checks

  if (!atom->angmom_flag || !atom->quat_flag ||
      !atom->avec->shape_type)
    error->all("Compute temp/asphere requires atom attributes "
	       "angmom, quat, shape");
  if (atom->radius_flag || atom->rmass_flag)
    error->all("Compute temp/asphere cannot be used with atom attributes "
	       "diameter or rmass");
}

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::~ComputeTempAsphere()
{
  delete [] id_bias;
  delete [] vector;
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (shape[type[i]][0] == 0.0)
	error->one("Compute temp/asphere requires extended particles");

  if (tempbias) {
    int i = modify->find_compute(id_bias);
    if (i < 0) error->all("Could not find compute ID for temperature bias");
    tbias = modify->compute[i];
    if (tbias->tempflag == 0)
      error->all("Bias compute does not calculate temperature");
    if (tbias->tempbias == 0)
      error->all("Bias compute does not calculate a velocity bias");
    if (tbias->igroup != igroup)
      error->all("Bias compute group does not match compute group");
    tbias->init();
    if (strcmp(tbias->style,"temp/region") == 0) tempbias = 2;
    else tempbias = 1;
  }

  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();

  calculate_inertia();
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::dof_compute()
{
  // 6 dof for 3d, 3 dof for 2d
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  double natoms = group->count(igroup);
  int nper = 6;
  if (domain->dimension == 2) nper = 3;
  dof = nper*natoms;

  // additional adjustments to dof

  if (tempbias == 1) dof -= tbias->dof_remove(-1) * natoms;
  else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
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

double ComputeTempAsphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **quat = atom->quat;
  double **angmom = atom->angmom;
  double *mass = atom->mass;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int itype;
  double wbody[3];
  double rot[3][3];
  double t = 0.0;
  
  // sum translationals and rotational energy for each particle
  // no point particles since divide by inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      itype = type[i];
      t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[itype];

      // wbody = angular velocity in body frame
      
      MathExtra::quat_to_mat(quat[i],rot);
      MathExtra::transpose_times_column3(rot,angmom[i],wbody);
      wbody[0] /= inertia[itype][0];
      wbody[1] /= inertia[itype][1];
      wbody[2] /= inertia[itype][2];

      t += inertia[itype][0]*wbody[0]*wbody[0] +
	inertia[itype][1]*wbody[1]*wbody[1] +
	inertia[itype][2]*wbody[2]*wbody[2];
    }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **quat = atom->quat;
  double **angmom = atom->angmom;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int itype;
  double wbody[3];
  double rot[3][3];
  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // translational kinetic energy

      itype = type[i];
      massone = mass[itype];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];

      // wbody = angular velocity in body frame

      MathExtra::quat_to_mat(quat[i],rot);
      MathExtra::transpose_times_column3(rot,angmom[i],wbody);
      wbody[0] /= inertia[itype][0];
      wbody[1] /= inertia[itype][1];
      wbody[2] /= inertia[itype][2];

      // rotational kinetic energy

      t[0] += inertia[itype][0]*wbody[0]*wbody[0];
      t[1] += inertia[itype][1]*wbody[1]*wbody[1];
      t[2] += inertia[itype][2]*wbody[2]*wbody[2];
      t[3] += inertia[itype][0]*wbody[0]*wbody[1];
      t[4] += inertia[itype][1]*wbody[0]*wbody[2];
      t[5] += inertia[itype][2]*wbody[1]*wbody[2];
    }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   principal moments of inertia for ellipsoids
------------------------------------------------------------------------- */

void ComputeTempAsphere::calculate_inertia()
{
  double *mass = atom->mass;
  double **shape = atom->shape;

  for (int i = 1; i <= atom->ntypes; i++) {
    inertia[i][0] = mass[i] * 
      (shape[i][1]*shape[i][1]+shape[i][2]*shape[i][2]) / 5.0;
    inertia[i][1] = mass[i] * 
      (shape[i][0]*shape[i][0]+shape[i][2]*shape[i][2]) / 5.0;
    inertia[i][2] = mass[i] * 
      (shape[i][0]*shape[i][0]+shape[i][1]*shape[i][1]) / 5.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempAsphere::remove_bias(int i, double *v)
{
  if (tbias) tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempAsphere::restore_bias(int i, double *v)
{
  if (tbias) tbias->restore_bias(i,v);
}
