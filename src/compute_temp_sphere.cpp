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

#include "mpi.h"
#include "string.h"
#include "compute_temp_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

#define INERTIA 0.4          // moment of inertia for sphere

/* ---------------------------------------------------------------------- */

ComputeTempSphere::ComputeTempSphere(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 3 && narg != 4)
    error->all("Illegal compute temp/sphere command");

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

  // error checks

  if (!atom->omega_flag)
    error->all("Compute temp/sphere requires atom attribute omega");
  if (!atom->radius_flag && !atom->avec->shape_type)
    error->all("Compute temp/sphere requires atom attribute "
	       "radius or shape");
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
  int i,itype;

  // if shape used, check that all particles are spherical
  // point particles are allowed

  if (atom->radius == NULL) {
    double **shape = atom->shape;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	itype = type[i];
	if (shape[itype][0] != shape[itype][1] || 
	    shape[itype][0] != shape[itype][2])
	  error->one("Compute temp/sphere requires "
		     "spherical particle shapes");
      }
  }

  if (tempbias) {
    i = modify->find_compute(id_bias);
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
  for (i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::dof_compute()
{
  int count,count_all;

  // 6 or 3 dof for extended/point particles for 3d
  // 3 or 2 dof for extended/point particles for 2d
  // assume full rotation of extended particles
  // user should correct this via compute_modify if needed

  int dimension = domain->dimension;

  double *radius = atom->radius;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  count = 0;
  if (dimension == 3) {
    if (radius) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (radius[i] == 0.0) count += 3;
	  else count += 6;
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (shape[type[i]][0] == 0.0) count += 3;
	  else count += 6;
	}
      }
    }
  } else {
    if (radius) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (radius[i] == 0.0) count += 2;
	  else count += 3;
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (shape[type[i]][0] == 0.0) count += 2;
	  else count += 3;
	}
      }
    }
  }

  MPI_Allreduce(&count,&count_all,1,MPI_INT,MPI_SUM,world);
  dof = count_all;

  // additional adjustments to dof

  if (tempbias == 1) {
    double natoms = group->count(igroup);
    dof -= tbias->dof_remove(-1) * natoms;

  } else if (tempbias == 2) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    count = 0;
    if (dimension == 3) {
      if (radius) {
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    if (tbias->dof_remove(i)) {
	      if (radius[i] == 0.0) count += 3;
	      else count += 6;
	    }
	  }
      } else {
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    if (tbias->dof_remove(i)) {
	      if (shape[type[i]][0] == 0.0) count += 3;
	      else count += 6;
	    }
	  }
      }
    } else {
      if (radius) {
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    if (tbias->dof_remove(i)) {
	      if (radius[i] == 0.0) count += 2;
	      else count += 3;
	    }
	  }
      } else {
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    if (tbias->dof_remove(i)) {
	      if (shape[type[i]][0] == 0.0) count += 2;
	      else count += 3;
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
  int i,itype;

  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // 4 cases depending on radius vs shape and rmass vs mass
  // point particles will not contribute rotation due to radius or shape = 0

  double t = 0.0;

  if (radius) {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	    rmass[i];
	  t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		omega[i][2]*omega[i][2]) * 
	    INERTIA*radius[i]*radius[i]*rmass[i];
	}

    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	    mass[itype];
	  t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		omega[i][2]*omega[i][2]) * 
	    INERTIA*radius[i]*radius[i]*mass[itype];
	}
    }

  } else {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	    rmass[i];
	  t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		omega[i][2]*omega[i][2]) *
	    INERTIA*shape[itype][0]*shape[itype][0]*rmass[i];
	}
      
    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * 
	    mass[itype];
	  t += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
		omega[i][2]*omega[i][2]) *
	    INERTIA*shape[itype][0]*shape[itype][0]*mass[itype];
	}
    }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic || tempbias == 2) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempSphere::compute_vector()
{
  int i,itype;

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
  }

  double **v = atom->v;
  double **omega = atom->omega;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // 4 cases depending on radius vs shape and rmass vs mass
  // point particles will not contribute rotation due to radius or shape = 0

  double massone,inertiaone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  if (radius) {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  massone = rmass[i];
	  t[0] += massone * v[i][0]*v[i][0];
	  t[1] += massone * v[i][1]*v[i][1];
	  t[2] += massone * v[i][2]*v[i][2];
	  t[3] += massone * v[i][0]*v[i][1];
	  t[4] += massone * v[i][0]*v[i][2];
	  t[5] += massone * v[i][1]*v[i][2];
	  
	  inertiaone = INERTIA*radius[i]*radius[i]*rmass[i];
	  t[0] += inertiaone * omega[i][0]*omega[i][0];
	  t[1] += inertiaone * omega[i][1]*omega[i][1];
	  t[2] += inertiaone * omega[i][2]*omega[i][2];
	  t[3] += inertiaone * omega[i][0]*omega[i][1];
	  t[4] += inertiaone * omega[i][0]*omega[i][2];
	  t[5] += inertiaone * omega[i][1]*omega[i][2];
	}

    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  massone = mass[itype];
	  t[0] += massone * v[i][0]*v[i][0];
	  t[1] += massone * v[i][1]*v[i][1];
	  t[2] += massone * v[i][2]*v[i][2];
	  t[3] += massone * v[i][0]*v[i][1];
	  t[4] += massone * v[i][0]*v[i][2];
	  t[5] += massone * v[i][1]*v[i][2];
	  
	  inertiaone = INERTIA*radius[i]*radius[i]*mass[itype];
	  t[0] += inertiaone * omega[i][0]*omega[i][0];
	  t[1] += inertiaone * omega[i][1]*omega[i][1];
	  t[2] += inertiaone * omega[i][2]*omega[i][2];
	  t[3] += inertiaone * omega[i][0]*omega[i][1];
	  t[4] += inertiaone * omega[i][0]*omega[i][2];
	  t[5] += inertiaone * omega[i][1]*omega[i][2];
	}
    }

  } else {
    if (rmass) {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  massone = rmass[i];
	  t[0] += massone * v[i][0]*v[i][0];
	  t[1] += massone * v[i][1]*v[i][1];
	  t[2] += massone * v[i][2]*v[i][2];
	  t[3] += massone * v[i][0]*v[i][1];
	  t[4] += massone * v[i][0]*v[i][2];
	  t[5] += massone * v[i][1]*v[i][2];
	  
	  inertiaone = INERTIA*shape[itype][0]*shape[itype][0]*rmass[i];
	  t[0] += inertiaone * omega[i][0]*omega[i][0];
	  t[1] += inertiaone * omega[i][1]*omega[i][1];
	  t[2] += inertiaone * omega[i][2]*omega[i][2];
	  t[3] += inertiaone * omega[i][0]*omega[i][1];
	  t[4] += inertiaone * omega[i][0]*omega[i][2];
	  t[5] += inertiaone * omega[i][1]*omega[i][2];
	}

    } else {
      for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  itype = type[i];
	  massone = mass[itype];
	  t[0] += massone * v[i][0]*v[i][0];
	  t[1] += massone * v[i][1]*v[i][1];
	  t[2] += massone * v[i][2]*v[i][2];
	  t[3] += massone * v[i][0]*v[i][1];
	  t[4] += massone * v[i][0]*v[i][2];
	  t[5] += massone * v[i][1]*v[i][2];
	  
	  inertiaone = INERTIA*shape[itype][0]*shape[itype][0]*mass[itype];
	  t[0] += inertiaone * omega[i][0]*omega[i][0];
	  t[1] += inertiaone * omega[i][1]*omega[i][1];
	  t[2] += inertiaone * omega[i][2]*omega[i][2];
	  t[3] += inertiaone * omega[i][0]*omega[i][1];
	  t[4] += inertiaone * omega[i][0]*omega[i][2];
	  t[5] += inertiaone * omega[i][1]*omega[i][2];
	}
    }
  }

  if (tempbias) tbias->restore_bias_all();

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempSphere::remove_bias(int i, double *v)
{
  tbias->remove_bias(i,v);
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempSphere::restore_bias(int i, double *v)
{
  tbias->restore_bias(i,v);
}
