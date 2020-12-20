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

FixBrownian::FixBrownian(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  virial_flag = 1;

  time_integrate = 1;

  if (narg != 6 && narg != 8)
    error->all(FLERR,"Illegal fix brownian command.");

  gamma_t = utils::numeric(FLERR,arg[3],false,lmp);
  if (gamma_t <= 0.0)
    error->all(FLERR,"Fix brownian viscous drag "
	       "coefficient must be > 0.");

  diff_t = utils::numeric(FLERR,arg[4],false,lmp);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix brownian diffusion "
	       "coefficient must be > 0.");
    
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
  if (seed <= 0) error->all(FLERR,"Fix brownian seed must be > 0.");

  noise_flag = 1;
  gaussian_noise_flag = 0;

  if (narg == 8) {
    
    if (strcmp(arg[6],"rng") == 0) {
      if (strcmp(arg[7],"uniform") == 0) {
	noise_flag = 1;
      } else if (strcmp(arg[7],"gaussian") == 0) {
	noise_flag = 1;
	gaussian_noise_flag = 1;
      } else if (strcmp(arg[7],"none") == 0) {
	noise_flag = 0;
      } else {
	error->all(FLERR,"Illegal fix brownian command.");
      }
    } else {
      error->all(FLERR,"Illegal fix brownian command.");
    }
  }
  
  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

int FixBrownian::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBrownian::~FixBrownian()
{
  delete random;
}



/* ---------------------------------------------------------------------- */

void FixBrownian::init()
{
  
  g1 =  force->ftm2v/gamma_t;
  if (noise_flag == 0) {
    g2 = 0;
    rng_func = &RanMars::zero_rng;
  } else if (gaussian_noise_flag == 1) {
    g2 = gamma_t*sqrt(2 * diff_t)/force->ftm2v;
    rng_func = &RanMars::gaussian;
  } else {
    g2 = gamma_t*sqrt( 24 * diff_t)/force->ftm2v;
    rng_func = &RanMars::uniform_middle;
  }

  dt = update->dt;
  sqrtdt = sqrt(dt);
}

void FixBrownian::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBrownian::initial_integrate(int /* vflag */)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


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
                  
    }
  }
  
  
  return;
}

/* ----------------------------------------------------------------------
   apply random force, stolen from MISC/fix_efield.cpp 
------------------------------------------------------------------------- */

void FixBrownian::post_force(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
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

void FixBrownian::reset_dt()
{

  dt = update->dt;
  sqrtdt = sqrt(dt);
}
