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
   Contributing authors: Carolyn Phillips (U Mich), reservoir energy tally
                         Aidan Thompson (SNL) GJF formulation
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_langevin_spin.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "math_const.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixLangevinSpin::FixLangevinSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_temp(NULL), random(NULL)
{
  if (narg != 7) error->all(FLERR,"Illegal fix langevin/spin command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  nevery = 1;

  temp = force->numeric(FLERR,arg[3]);
  alpha_t = force->numeric(FLERR,arg[4]);
  alpha_l = force->numeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);

  if (alpha_t < 0.0) error->all(FLERR,"Fix langevin/spin transverse damping must be >= 0.0");
  if (alpha_l < 0.0) error->all(FLERR,"Fix langevin/spin transverse damping must be >= 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix langevin/spin seed must be > 0");

  // initialize Marsaglia RNG with processor-unique seed
  //random = new RanMars(lmp,seed + comm->me);
  random = new RanPark(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixLangevinSpin::~FixLangevinSpin()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixLangevinSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::init()
{
  // warn if any fix comes after this one  
  int after = 0;
  int flag_force = 0;
  int flag_lang = 0;
  for (int i = 0; i < modify->nfix; i++) { 
     if (strcmp("force/spin",modify->fix[i]->style)==0) flag_force = MAX(flag_force,i);
     if (strcmp("langevin/spin",modify->fix[i]->style)==0) flag_lang = i;
  }
  if (flag_force >= flag_lang) error->all(FLERR,"Fix langevin/spin should come after all other spin fixes");  

  dts = update->dt; 
  Gil_factor = 1.0/(1.0+(alpha_t)*(alpha_t));
  
  double hbar = force->hplanck/MY_2PI; //eV/(rad.THz)
  double kb = force->boltz;
  D = (MY_2PI*Gil_factor*kb*temp)/hbar/dts;
  sigma = sqrt(D);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::post_force(int vflag)
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;  
              
  double sx, sy, sz;
  double fmx, fmy, fmz;
  double cpx, cpy, cpz;
  double rx, ry, rz;  
                          
  // apply transverse magnetic damping to spins
  // add the damping to the effective field
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
		sx = sp[i][0];//Getting Mag. and Mag. force components
		sy = sp[i][1];
		sz = sp[i][2];
		
		fmx = fm[i][0];
		fmy = fm[i][1];
		fmz = fm[i][2];
	
		cpx = fmy*sz - fmz*sy;//Computing cross product
		cpy = fmz*sx - fmx*sz;
		cpz = fmx*sy - fmy*sx;
		
		fmx -= alpha_t*cpx;//Taking the damping value away
		fmy -= alpha_t*cpy;
		fmz -= alpha_t*cpz;
		
		fm[i][0] = fmx;
		fm[i][1] = fmy;
		fm[i][2] = fmz;		
   }

  //printf("test damping. 1;i=0, fx=%g, fy=%g, fz=%g \n",fm[0][0],fm[0][1],fm[0][2]);
  //apply thermal effects
  //add random field to fm 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
		rx = sigma*random->gaussian();//Drawing random distributions
		ry = sigma*random->gaussian();
		rz = sigma*random->gaussian();
	
		//rx = sigma*(random->uniform() - 0.5);
		//ry = sigma*(random->uniform() - 0.5);
		//rz = sigma*(random->uniform() - 0.5);
	
        	fm[i][0] += rx;//Adding random field
		fm[i][1] += ry;
		fm[i][2] += rz;
                
                fm[i][0] *= Gil_factor;//Multiplying by Gilbert's prefactor 
                fm[i][1] *= Gil_factor; 
                fm[i][2] *= Gil_factor; 
		
   }

   //printf("test langevin 1;i=0, fx=%g, fy=%g, fz=%g \n",fm[0][0],fm[0][1],fm[0][2]);

   //printf("test rand var: %g, sigma=%g \n",(random->uniform()-0.5),sigma);
   //printf("test rand var: %g, sigma=%g \n",random->gaussian(),sigma);
   //printf("test random 1;i=0, fx=%g, fy=%g, fz=%g \n",fm[0][0],fm[0][1],fm[0][2]);  
   //printf("test dt: %g, sigma=%g \n",dts,random->gaussian(),sigma);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

