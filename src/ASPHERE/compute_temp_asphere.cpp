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
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::ComputeTempAsphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute temp command");

  if (!atom->quat_flag || !atom->angmom_flag)
    error->all("Compute temp/asphere requires atom attributes quat, angmom");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_temp_sphere:inertia");
}

/* ---------------------------------------------------------------------- */

ComputeTempAsphere::~ComputeTempAsphere()
{
  delete [] vector;
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  recount();

  calculate_inertia();
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::recount()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;

  // add rotational degrees of freedom
  // 0 for sphere, 2 for uniaxial, 3 for biaxial

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  int itype;
  int rot_dof = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      itype = type[i];
      if (shape[itype][0] == shape[itype][1] && 
	  shape[itype][1] == shape[itype][2]) continue;
      else if (shape[itype][0] == shape[itype][1] || 
	       shape[itype][1] == shape[itype][2] ||
	       shape[itype][0] == shape[itype][2]) rot_dof += 2;
      else rot_dof += 3;
    }

  int rot_total;
  MPI_Allreduce(&rot_dof,&rot_total,1,MPI_INT,MPI_SUM,world);
  dof += rot_total;
    
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempAsphere::compute_scalar()
{
  invoked |= INVOKED_SCALAR;

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
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // translational kinetic energy
      
      itype = type[i];
      t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[itype];

      // wbody = angular velocity in body frame

      if (!(shape[itype][0] == shape[itype][1] && 
            shape[itype][1] == shape[itype][2])) {

        MathExtra::quat_to_mat(quat[i],rot);
        MathExtra::transpose_times_column3(rot,angmom[i],wbody);
        wbody[0] /= inertia[itype][0];
        wbody[1] /= inertia[itype][1];
        wbody[2] /= inertia[itype][2];

        // rotational kinetic energy

        t += inertia[itype][0]*wbody[0]*wbody[0]+
             inertia[itype][1]*wbody[1]*wbody[1]+
             inertia[itype][2]*wbody[2]*wbody[2];
      }
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) recount();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempAsphere::compute_vector()
{
  int i;

  invoked |= INVOKED_VECTOR;

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
