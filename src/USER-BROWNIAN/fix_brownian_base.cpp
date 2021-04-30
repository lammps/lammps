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

#include "fix_brownian.h"

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

/* ---------------------------------------------------------------------- */

FixBrownianBase::FixBrownianBase(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  time_integrate = 1;

  noise_flag = 1;
  gaussian_noise_flag = 0;
  gamma_t_flag = gamma_r_flag = 0;
  gamma_t_eigen_flag = gamma_r_eigen_flag = 0;
  dipole_flag = 0;
  
  if (narg < 5)
    error->all(FLERR,"Illegal fix brownian command.");

  temp = utils::numeric(FLERR,arg[3],false,lmp);
  if (temp <= 0) error->all(FLERR,"Fix brownian temp must be > 0.");
  
  seed = utils::inumeric(FLERR,arg[4],false,lmp);
  if (seed <= 0) error->all(FLERR,"Fix brownian seed must be > 0.");


  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"rng") == 0) {
      if (narg == iarg + 1) {
	error->all(FLERR,"Illegal fix brownian command.");
      }
      if (strcmp(arg[iarg + 1],"uniform") == 0) {
	noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"gaussian") == 0) {
	noise_flag = 1;
	gaussian_noise_flag = 1;
      } else if (strcmp(arg[iarg + 1],"none") == 0) {
	noise_flag = 0;
      } else {
	error->all(FLERR,"Illegal fix brownian command.");
      }
      iarg = iarg + 2;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      if (narg == iarg + 3) {
	error->all(FLERR,"Illegal fix brownian command.");
      }

      dipole_flag = 1;
      dipole_body = new double[3];
      
      dipole_body[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      dipole_body[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      dipole_body[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg = iarg + 4;
      
    } else if (strcmp(arg[iarg],"gamma_t_eigen") == 0) {
      if (narg == iarg + 3) {
	error->all(FLERR,"Illegal fix brownian command.");
      }

      gamma_t_eigen_flag = 1;
      gamma_t_inv = new double[3];
      gamma_t_invsqrt = new double[3];
      gamma_t_inv[0] = 1./utils::numeric(FLERR,arg[iarg+1],false,lmp);
      gamma_t_inv[1] = 1./utils::numeric(FLERR,arg[iarg+2],false,lmp);

      if (domain->dimension == 2) {
	if (strcmp(arg[iarg+3],"inf") != 0) {
	  error->all(FLERR,"Fix brownian gamma_t_eigen third value must be inf for 2D system.");
	}
	gamma_t_inv[2] = 0;
      } else { 
	gamma_t_inv[2] = 1./utils::numeric(FLERR,arg[iarg+3],false,lmp);
      }

      if (gamma_t_inv[0] < 0 || gamma_t_inv[1] < 0 || gamma_t_inv[2] < 0) 
	error->all(FLERR,"Fix brownian gamma_t_eigen values must be > 0.");      
      
      gamma_t_invsqrt[0] = sqrt(gamma_t_inv[0]);
      gamma_t_invsqrt[1] = sqrt(gamma_t_inv[1]);
      gamma_t_invsqrt[2] = sqrt(gamma_t_inv[2]);      
      iarg = iarg + 4;

    } else if (strcmp(arg[iarg],"gamma_r_eigen") == 0) {
      if (narg == iarg + 3) {
	error->all(FLERR,"Illegal fix brownian command.");
      }

      gamma_r_eigen_flag = 1;      
      gamma_r_inv = new double[3];
      gamma_r_invsqrt = new double[3];


      if (domain->dimension == 2) {
	if (strcmp(arg[iarg+1],"inf") != 0) {
	  error->all(FLERR,"Fix brownian gamma_r_eigen first value must be inf for 2D system.");
	}
	gamma_r_inv[0] = 0;

	if (strcmp(arg[iarg+2],"inf") != 0) {
	  error->all(FLERR,"Fix brownian gamma_r_eigen second value must be inf for 2D system.");
	}
	gamma_r_inv[1] = 0;
      } else {

	gamma_r_inv[0] = 1./utils::numeric(FLERR,arg[iarg+1],false,lmp);
	gamma_r_inv[1] = 1./utils::numeric(FLERR,arg[iarg+2],false,lmp);

      }
      
      gamma_r_inv[2] = 1./utils::numeric(FLERR,arg[iarg+3],false,lmp);

      if (gamma_r_inv[0] < 0 || gamma_r_inv[1] < 0 || gamma_r_inv[2] < 0) 
	error->all(FLERR,"Fix brownian gamma_r_eigen values must be > 0.");
      
      gamma_r_invsqrt[0] = sqrt(gamma_r_inv[0]);
      gamma_r_invsqrt[1] = sqrt(gamma_r_inv[1]);
      gamma_r_invsqrt[2] = sqrt(gamma_r_inv[2]);
      iarg = iarg + 4;
      
    } else if (strcmp(arg[iarg],"gamma_t") == 0) {
      if (narg == iarg + 1) {
	error->all(FLERR,"Illegal fix brownian command.");
      }
      
      gamma_t_flag = 1;
      gamma_t = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (gamma_t <= 0) 
	error->all(FLERR,"Fix brownian gamma_t must be > 0.");
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg],"gamma_r") == 0) {
      if (narg == iarg + 1) {
	error->all(FLERR,"Illegal fix brownian command.");
      }

      gamma_r_flag = 1;
      gamma_r = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (gamma_r <= 0) 
	error->all(FLERR,"Fix brownian gamma_r must be > 0.");
      iarg = iarg + 2;

    } else {
      error->all(FLERR,"Illegal fix brownian command.");
    }
  }
  
  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

int FixBrownianBase::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBrownianBase::~FixBrownianBase()
{

  if (gamma_t_eigen_flag) {
    delete [] gamma_t_inv;
    delete [] gamma_t_invsqrt;
  }
  if (gamma_r_eigen_flag) {
    delete [] gamma_r_inv;
    delete [] gamma_r_invsqrt;
  }

  if (dipole_flag) {
    delete [] dipole_body;
  }

  
  delete random;
}



/* ---------------------------------------------------------------------- */

void FixBrownianBase::init()
{
  dt = update->dt;
  sqrtdt = sqrt(dt);
  
  g1 =  force->ftm2v;
  if (noise_flag == 0) {
    g2 = 0;
  } else if (gaussian_noise_flag == 1) {
    g2 = sqrt(2 * force->boltz*temp/dt/force->mvv2e);
  } else {
    g2 = sqrt( 24 * force->boltz*temp/dt/force->mvv2e);
  }


}



void FixBrownianBase::reset_dt()
{
  double sqrtdt_old = sqrtdt;
  dt = update->dt;
  sqrtdt = sqrt(dt);
  g2 *= sqrtdt_old/sqrtdt;
  
}
