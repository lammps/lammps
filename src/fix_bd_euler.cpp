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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_bd_euler.h"
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBDEuler::FixBDEuler(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix bd/euler command");

  if (domain->dimension != 2) error->all(FLERR,"Fix bd/euler requires 2d simulation");

  t_start = force->numeric(FLERR,arg[3]);
  t_target = t_start;
  t_stop = force->numeric(FLERR,arg[4]);
  diff = force->numeric(FLERR,arg[5]);
  if (diff <= 0.0) error->all(FLERR,"Fix bd/euler diffusion coefficient must be > 0.0");
  seed = force->inumeric(FLERR,arg[6]);
  if (seed <= 0) error->all(FLERR,"Illegal fix bd/euler command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixBDEuler::~FixBDEuler()
{

  delete random;

}


/* ---------------------------------------------------------------------- */

void FixBDEuler::init()
{

  if (!atom->mu_flag)
    error->all(FLERR,"Fix bd/euler requires atom attributes mu");

  //TBD: check if muz is zero for all particles


  // set square root of temperature
  compute_target();

  gamma1 = diff / force->boltz; 
  gamma2 = sqrt( 24 * diff );
  gamma3 = sqrt( 24 * 3 * diff );
  gamma4 = 3 * diff / force->boltz; 

  FixNVE::init();
}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixBDEuler::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}


/* ---------------------------------------------------------------------- */

void FixBDEuler::initial_integrate(int vflag)
{
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **v = atom->v;
  double **mu = atom->mu;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **torque = atom->torque;
  double **omega = atom->omega;


  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);

  // set square root of temperature
  compute_target();

  //printf ("\n\n\n");
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {      

      //printf("%d: (%.18f %.8f %.8f) (%.8f %.8f %.8f) \n",i , x[i][0], x[i][1], x[i][2],  mu[i][0], mu[i][1], mu[i][2] );
      

      da = dt * gamma1 * f[i][0] / t_target   +    sqrtdt * gamma2 * (random->uniform()-0.5);
      double dx = da;
      x[i][0] +=  da;
      v[i][0]  =  da/dt;
      da = dt * gamma1 * f[i][1] / t_target   +    sqrtdt * gamma2 * (random->uniform()-0.5);
      x[i][1] +=  da;
      double dy = da;
      v[i][1]  =  da/dt;
      //printf("   (%.18f %.8f %.8f) (%.8f %.8f %.8f) \n", x[i][0], x[i][1], x[i][2],  mu[i][0], mu[i][1], mu[i][2] );

      dar   =   sqrtdt * gamma3 * (random->uniform()-0.5);
      da    =   dt * gamma4 * torque[i][2] / t_target + dar;
      cosda = cos(da);
      sinda = sin(da);
      omega[i][2] = da/dt; 


      //printf("%d: x=(%f %f %f) mu=(%f %f %f) torque=(%f %f %f) f=(%f %f %f) \n",i , x[i][0], x[i][1], x[i][2],  mu[i][0], mu[i][1], mu[i][2] , torque[i][0], torque[i][1], torque[i][2], f[i][0], f[i][1], f[i][2]);
      //printf("      da_r=%f, da_torque=%f, d_theta=%f, theta=%f\n\n", dar/3.1415*180, (da-dar)/3.1415*180, da/3.1415*180, atan2(mu[i][1], mu[i][0])/3.1415*180 );

      //printf("d_theta=%.8f, cos=%.8f, sin=%.8f\n", da/3.1415*180, cosda, sinda );

      //printf("%d: dx=(%g %g) v=(%g %g) dtheta=%g  omega=%g \n",i, dx,dy, v[i][0],v[i][1], da, omega[i][2] );

      da = mu[i][0];
      mu[i][0] =  mu[i][0]*cosda - mu[i][1]*sinda;
      mu[i][1] =  sinda * da     + mu[i][1]*cosda;
      //printf("   (%.18f %.8f %.8f) (%.8f %.8f %.8f) \n", x[i][0], x[i][1], x[i][2],  mu[i][0], mu[i][1], mu[i][2] );
      //printf("   dipole_squered = %.8f \n\n", mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1]);

      //torque[i][2] = 0.0;
      //torque[i][1] = 0.0;




    }

}


