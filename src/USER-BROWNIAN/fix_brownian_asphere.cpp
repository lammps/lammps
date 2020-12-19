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

#include "fix_brownian_asphere.h"

#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
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

enum{NODIPOLE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixBrownianAsphere::FixBrownianAsphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  virial_flag = 1;

  time_integrate = 1;

  dipole_flag = NODIPOLE;


  if (narg > 11 || narg < 8 )
    error->all(FLERR,"Illegal fix brownian/asphere command.");

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix brownian/asphere requires atom style sphere");

  gamma_t = utils::numeric(FLERR,arg[3],false,lmp);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix brownian/asphere translational viscous drag "
	       "coefficient must be > 0.");

  gamma_r = utils::numeric(FLERR,arg[4],false,lmp);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix brownian/asphere rotational viscous drag "
	       "coefficient must be > 0.");


  diff_t = utils::numeric(FLERR,arg[5],false,lmp);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix brownian/asphere translational diffusion "
	       "coefficient must be > 0.");
  
  diff_r = utils::numeric(FLERR,arg[6],false,lmp);
  if (diff_r <= 0.0)
    error->all(FLERR,"Fix brownian/asphere rotational diffusion "
	       "coefficient must be > 0.");
  
  seed = utils::inumeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Fix brownian/asphere seed must be > 0.");

  noise_flag = 1;
  gaussian_noise_flag = 0;

  int iarg = 8;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rng") == 0) {
      if (narg == iarg + 1) {
	error->all(FLERR,"Illegal fix/brownian/asphere command.");
      }
      if (strcmp(arg[iarg + 1],"uniform") == 0) {
	noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"gaussian") == 0) {
	noise_flag = 1;
	gaussian_noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"none") == 0) {
	noise_flag = 0;
      } else {
	error->all(FLERR,"Illegal fix/brownian/asphere command.");
      }
      iarg = iarg + 2;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      dipole_flag = DIPOLE;
      iarg = iarg + 1;
    } else {
      error->all(FLERR,"Illegal fix/brownian/asphere command.");
    }
  }
  
  if (dipole_flag == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix brownian/asphere dipole requires atom attribute mu");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

int FixBrownianAsphere::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBrownianAsphere::~FixBrownianAsphere()
{
  delete random;
}



/* ---------------------------------------------------------------------- */

void FixBrownianAsphere::init()
{


  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute brownian/asphere requires "
	       "atom style ellipsoid");
  
  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK
  
  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
	error->one(FLERR,"Fix brownian/asphere requires extended particles");
  
  
  if (dipole_flag == DIPOLE) {
    
    double f_act[3] = { 1.0, 0.0, 0.0 };
    double f_rot[3];
    double *quat;
    int *ellipsoid = atom->ellipsoid;
    AtomVecEllipsoid::Bonus *bonus = avec->bonus;
    
    double Q[3][3];
    
    double **mu = atom->mu;
    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	quat = bonus[ellipsoid[i]].quat;
	MathExtra::quat_to_mat( quat, Q );
	MathExtra::matvec( Q, f_act, f_rot );
	
	mu[i][0] = f_rot[0];
	mu[i][1] = f_rot[1];
	mu[i][2] = f_rot[2];
	
      }
    }   
  }

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

void FixBrownianAsphere::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */


void FixBrownianAsphere::initial_integrate(int /* vflag */)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  int d3rot;  // whether to compute angular momentum in xy plane
  
  if (domain->dimension==2) {
    d3rot = 0;
  } else {
    d3rot = 1;
  }

  if (dipole_flag == DIPOLE) {

    // if dipole is being tracked, then update it along with
    // quaternions accordingly along with angular velocity

    double wq[4];
    
    double *quat;
    int *ellipsoid = atom->ellipsoid;
    AtomVecEllipsoid::Bonus *bonus = avec->bonus;
    
    // project dipole along x axis of quat
    double f_act[3] = { 1.0, 0.0, 0.0 }; 
    double f_rot[3];

    double Q[3][3];
    
    double **mu = atom->mu;
    

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	
	update_x_and_omega(x[i],v[i],omega[i],f[i],torque[i],d3rot);

	quat = bonus[ellipsoid[i]].quat;
	
	MathExtra::vecquat(omega[i],quat,wq);
	
	quat[0] = quat[0] + 0.5*dt*wq[0];
	quat[1] = quat[1] + 0.5*dt*wq[1];
	quat[2] = quat[2] + 0.5*dt*wq[2];
	quat[3] = quat[3] + 0.5*dt*wq[3];
	MathExtra::qnormalize(quat);

	MathExtra::quat_to_mat( quat, Q );
	MathExtra::matvec( Q, f_act, f_rot );

	mu[i][0] = f_rot[0];
	mu[i][1] = f_rot[1];
	mu[i][2] = f_rot[2];

      }
    }
  } else {

    // if no dipole, just update quaternions and
    // angular velocity


    double wq[4];

    double *quat;
    int *ellipsoid = atom->ellipsoid;
    AtomVecEllipsoid::Bonus *bonus = avec->bonus;
    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {

	update_x_and_omega(x[i],v[i],omega[i],f[i],torque[i],d3rot);
	
	quat = bonus[ellipsoid[i]].quat;
	
	MathExtra::vecquat(omega[i],quat,wq);
	
	quat[0] = quat[0] + 0.5*dt*wq[0];
	quat[1] = quat[1] + 0.5*dt*wq[1];
	quat[2] = quat[2] + 0.5*dt*wq[2];
	quat[3] = quat[3] + 0.5*dt*wq[3];
	MathExtra::qnormalize(quat);

      }
    }
  }
  return;
}

void FixBrownianAsphere::update_x_and_omega(double *x, double *v, double *omega,
				   double *f, double *torque, int d3rot)
{
  double dx, dy, dz;
  
  dx = dt * g1 * f[0];
  x[0] +=  dx;
  v[0]  =  dx/dt;
  
  dy = dt * g1 * f[1];
  x[1] +=  dy;
  v[1]  =  dy/dt;
  
  dz = dt * g1 * f[2];
  x[2] +=  dz;
  v[2]  =  dz/dt;
  
  omega[0] = d3rot * g3* torque[0];
  omega[1] = d3rot * g3* torque[1];
  omega[2] = g3* torque[2];

  return;
}

/* ----------------------------------------------------------------------
   apply random force, stolen from MISC/fix_efield.cpp 
------------------------------------------------------------------------- */

void FixBrownianAsphere::post_force(int vflag)
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

void FixBrownianAsphere::reset_dt()
{

  dt = update->dt;
  sqrtdt = sqrt(dt);
}
