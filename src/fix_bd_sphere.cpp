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
   Originally modified from USER-CGDNA/fix_nve_dotc_langevin.cpp. 

   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_bd_sphere.h"
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

FixBdSphere::FixBdSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  time_integrate = 1;
  
  if (narg != 8 && narg != 10 && narg != 12)
    error->all(FLERR,"Illegal fix bd/sphere command.");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix bd/sphere requires atom style sphere");
  if (!atom->mu_flag)
    error->all(FLERR,"Fix bd/sphere requires atom attribute mu");

  gamma_t = force->numeric(FLERR,arg[3]);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere translational viscous drag "
	       "coefficient must be > 0.");

  gamma_r = force->numeric(FLERR,arg[4]);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere rotational viscous drag "
	       "coefficient must be > 0.");

  diff_t = force->numeric(FLERR,arg[5]);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere translational diffusion "
	       "coefficient must be > 0.");
  
  diff_r = force->numeric(FLERR,arg[6]);
  if (diff_r <= 0.0)
    error->all(FLERR,"Fix bd/sphere rotational diffusion "
	       "coefficient must be > 0.");
  
  seed = force->inumeric(FLERR,arg[7]);
  if (seed <= 0) error->all(FLERR,"Fix bd/sphere seed must be > 0.");

  noise_flag = 1;
  gaussian_noise_flag = 0;
  rotate_planar_flag = 0;

  int iarg == 8;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rng") == 0) {
      if (strcmp(arg[iarg + 1],"uniform") == 0) {
	noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"gaussian") == 0) {
	noise_flag = 1;
	gaussian_noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"none") == 0) {
	noise_flag = 0;
      } else {
	error->all(FLERR,"Illegal fix/bd/sphere command.");
      }
    } else if (strcmp(arg[iarg],"rotate_planar") == 0) {

      if (strcmp(arg[iarg + 1],"yes") == 0) {
	rotate_planar_flag = 1;
	if (domain->dimension != 2) {
	  error->all(FLERR,"Cannot constrain rotational degrees of freedom "
		     "to the xy plane if the simulation is in 3D "
		     "(in fix/bd/sphere).");
	}
      } else if (strcmp(arg[iarg + 1],"no") == 0) {
	rotate_planar_flag = 0;
      } else {
	error->all(FLERR,"Illegal fix/bd/sphere command.");
      }
    } else {
      error->all(FLERR,"Illegal fix/bd/sphere command.");
    }
    iarg = iarg + 2;
  }

  
  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);


}

/* ---------------------------------------------------------------------- */

int FixBdSphere::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBdSphere::~FixBdSphere()
{

  delete random;

}


/* ---------------------------------------------------------------------- */

void FixBdSphere::init()
{

  g1 =  force->ftm2v/gamma_t;
  g3 = force->ftm2v/gamma_r;
  if (noise_flag == 0) {
    g2 = 0;
    g4 = 0;
    rng_func = &RanMars::zero_rng;
  } else if (gaussian_noise_flag == 1) {
    g2 = sqrt(2 * diff_t);
    g4 = sqrt(2 * diff_r);
    rng_func = &RanMars::gaussian;
  } else {
    g2 = sqrt( 24 * diff_t);
    g4 = sqrt( 24 * diff_r );
    rng_func = &RanMars::uniform_middle;
  }

  if (domain->dimension == 2 && rotate_planar_flag == 0) {
    error->warning(FLERR,"Using a 2D simulation, but allowing for "
		   "full (3D) rotation (in fix/bd/sphere).");
  }

  
}

/* ---------------------------------------------------------------------- */

void FixBdSphere::initial_integrate(int /* vflag */)
{
  double **x = atom->x;
  double **v = atom->v;
  double **mu = atom->mu;
  double **f = atom->f;
  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  double dtheta;
  double mux,muy,muz,mu_tmp,wx,wy,wz;
  double prefac_1, prefac_2;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed 
  dt = update->dt;
  sqrtdt = sqrt(dt);

  int d3rot;  // whether to compute angular momentum in xy plane
  
  if (rotate_planar_flag) {
    d3rot = 0;
  } else {
    d3rot = 1;
  }

  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      
      dx = (dt * g1 * f[i][0]
	    +    sqrtdt * g2 * (random->*rng_func)());
      x[i][0] +=  dx;
      v[i][0]  =  dx/dt;
      dy = (dt * g1 * f[i][1] 
	    +    sqrtdt * g2 * (random->*rng_func)());
      x[i][1] +=  dy;
      v[i][1]  =  dy/dt;
      
      dz = (dt * g1 * f[i][2] 
	    +    sqrtdt * g2 * (random->*rng_func)());
      x[i][2] +=  dz;
      v[i][2]  =  dz/dt;
      
      
      omega[i][0] = d3rot*(g3* torque[i][0]
			   +  g4 * (random->*rng_func)()/sqrtdt);
      omega[i][1] = d3rot*(g3* torque[i][1]
			   +  g4 * (random->*rng_func)()/sqrtdt);
      omega[i][2] = (g3* torque[i][2]
		     + g4 * (random->*rng_func)()/sqrtdt);
      
      dtheta = sqrt((omega[i][0]*dt)**2+(omega[i][1]*dt)**2+(omega[i][2]*dt)**2);
      
      if (abs(dtheta) < 1e-14) {
	prefac_1 = dt;
	prefac_2 = 0.5*dt*dt;
      } else {
	prefac_1 = dt*sin(dtheta)/dtheta;
	prefac_2 = dt*dt*(1-cos(dtheta))/(dtheta*dtheta);
      }
      
      mux = mu[i][0];
      muy = mu[i][1];
      muz = mu[i][2];
      
      wx = omega[i][0];
      wy = omega[i][1];
      wz = omega[i][2];
      
      mu[i][0] = (mux + prefac_1 * ( -wz*muy + wy*muz )
		  + prefac_2 * ( -1*( wz*wz + wy*wy ) * mux
				 + ( wz*muz + wy*muy ) * wx));
      
      mu[i][1] = (muy + prefac_1 * ( wz*mux - wx*muz )
		  + prefac_2 * ( -1*(wz*wz + wx*wx) * muy
				 + ( wz*muz + wx*mux ) * wy));
      
      mu[i][2] = (muz + prefac_1 * ( -wy*mux + wx*muy )
		  + prefac_2 * ( -1*( wx*wx + wy*wy ) * muz
				 + ( wy*muy + wx*mux ) * wz));

      mu_tmp = sqrt(mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2]);

      mu[i][0] = mu[i][0]/mu_tmp;
      mu[i][1] = mu[i][1]/mu_tmp;
      mu[i][2] = mu[i][2]/mu_tmp;
      
    }
  }
  
  return;
}
