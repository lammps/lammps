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

#include "fix_bd_sphere.h"

#include <cmath>
#include <cstring>
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

#define SMALL 1e-14

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixBdSphere::FixBdSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  virial_flag = 1;

  time_integrate = 1;

  extra = NONE;

  if (narg > 11 || narg < 8 )
    error->all(FLERR,"Illegal fix bd/sphere command.");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix bd/sphere requires atom style sphere");

  gamma_t = utils::numeric(FLERR,arg[3],false,lmp);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere translational viscous drag "
	       "coefficient must be > 0.");

  gamma_r = utils::numeric(FLERR,arg[4],false,lmp);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere rotational viscous drag "
	       "coefficient must be > 0.");


  diff_t = utils::numeric(FLERR,arg[5],false,lmp);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix bd/sphere translational diffusion "
	       "coefficient must be > 0.");
  
  diff_r = utils::numeric(FLERR,arg[6],false,lmp);
  if (diff_r <= 0.0)
    error->all(FLERR,"Fix bd/sphere rotational diffusion "
	       "coefficient must be > 0.");
  
  seed = utils::inumeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Fix bd/sphere seed must be > 0.");

  noise_flag = 1;
  gaussian_noise_flag = 0;

  int iarg = 8;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rng") == 0) {
      if (narg == iarg + 1) {
	error->all(FLERR,"Illegal fix/bd/sphere command.");
      }
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
      iarg = iarg + 2;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      extra = DIPOLE;
      iarg = iarg + 1;
    } else {
      error->all(FLERR,"Illegal fix/bd/sphere command.");
    }
  }
  
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix bd/sphere update dipole requires atom attribute mu");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

int FixBdSphere::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
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
    g2 = gamma_t*sqrt(2 * diff_t)/force->ftm2v;
    g4 = gamma_r*sqrt(2 * diff_r)/force->ftm2v;
    rng_func = &RanMars::gaussian;
  } else {
    g2 = gamma_t*sqrt( 24 * diff_t)/force->ftm2v;
    g4 = gamma_r*sqrt( 24 * diff_r )/force->ftm2v;
    rng_func = &RanMars::uniform_middle;
  }

  dt = update->dt;
  sqrtdt = sqrt(dt);
}

void FixBdSphere::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBdSphere::initial_integrate(int /* vflag */)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  int d3rot;  // whether to compute angular momentum in xy plane

  dt = update->dt;
  sqrtdt = sqrt(dt);
  
  if (domain->dimension==2) {
    d3rot = 0;
  } else {
    d3rot = 1;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      
      dx = dt * g1 * f[i][0];
      x[i][0] +=  dx;
      v[i][0]  =  dx/dt;
      
      dy = dt * g1 * f[i][1];
      x[i][1] +=  dy;
      v[i][1]  =  dy/dt;
      
      dz = dt * g1 * f[i][2];
      x[i][2] +=  dz;
      v[i][2]  =  dz/dt;
      
      omega[i][0] = d3rot * g3* torque[i][0];
      omega[i][1] = d3rot * g3* torque[i][1];
      omega[i][2] = g3* torque[i][2];
            
    }
  }
  
  if (extra == DIPOLE) {
    
    double **mu = atom->mu;
    double dtheta;
    double mux,muy,muz,mu_tmp,wx,wy,wz;
    double prefac_1, prefac_2;
    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	
	dtheta = sqrt((omega[i][0]*dt)*(omega[i][0]*dt)
		      +(omega[i][1]*dt)*(omega[i][1]*dt)
		      +(omega[i][2]*dt)*(omega[i][2]*dt));
	
	
	if (fabs(dtheta) < SMALL) {
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
  }
  
  return;
}

/* ----------------------------------------------------------------------
   apply random force, stolen from MISC/fix_efield.cpp 
------------------------------------------------------------------------- */

void FixBdSphere::post_force(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  double **torque = atom->torque;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  
  // virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  double fx,fy,fz;
  double v[6];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      
      fx = g2 * (random->*rng_func)()/sqrtdt;
      fy = g2 * (random->*rng_func)()/sqrtdt;
      fz = g2 * (random->*rng_func)()/sqrtdt;
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;

      torque[i][0] = g4*(random->*rng_func)()/sqrtdt;
      torque[i][1] = g4*(random->*rng_func)()/sqrtdt;
      torque[i][2] = g4*(random->*rng_func)()/sqrtdt;

	if (evflag) {
	  v[0] = fx*x[i][0];
	  v[1] = fy*x[i][1];
	  v[2] = fz*x[i][2];
	  v[3] = fx*x[i][1];
	  v[4] = fx*x[i][2];
	  v[5] = fy*x[i][2];
	  v_tally(i, v);
	}
    }
}
