/* -------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
    Contributing authors: 
             Rodrigo Freitas   (Unicamp/Brazil) - rodrigohb@gmail.com
             Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
------------------------------------------------------------------------- */
 
#include "stdlib.h" 
#include "string.h"  
#include "math.h"
#include "fix_ti_rs.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

// Class constructor initialize all variables.

FixTIRS::FixTIRS(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) 
{
  // Checking the input information.
  if (narg < 7 || narg > 9) error->all(FLERR,"Illegal fix ti/rs command");

  // Fix flags.
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector   = 1;

  // Time variables.
  t_switch  = atoi(arg[5]);
  t_equil   = atoi(arg[6]);
  t0 = update->ntimestep;    
  if (t_switch < 0.0) error->all(FLERR,"Illegal fix ti/rs command");
  if (t_equil  < 0.0) error->all(FLERR,"Illegal fix ti/rs command");
 
  // Coupling parameter limits and initialization.
  l_initial = atof(arg[3]);
  l_final   = atof(arg[4]);
  sf = 1;
  if (narg > 7) {
    if (strcmp(arg[7], "function") == 0) sf = atoi(arg[8]);
    else error->all(FLERR,"Illegal fix ti/rs switching function");
    if ((sf<1) || (sf>3))
      error->all(FLERR,"Illegal fix ti/rs switching function");
  }
  lambda  =  switch_func(0);
  dlambda = dswitch_func(0);
}

/* ---------------------------------------------------------------------- */

FixTIRS::~FixTIRS() 
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);
}

/* ---------------------------------------------------------------------- */

int FixTIRS::setmask() 
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTIRS::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTIRS::setup(int vflag) 
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

void FixTIRS::min_setup(int vflag)
{
  post_force(vflag);
}


/* ---------------------------------------------------------------------- */

void FixTIRS::post_force(int vflag) 
{

  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  double **f = atom->f;

  // Scaling forces.
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] = lambda * f[i][0]; 
      f[i][1] = lambda * f[i][1];
      f[i][2] = lambda * f[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTIRS::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTIRS::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTIRS::initial_integrate(int vflag) 
{
  // Update the coupling parameter value.
  double t = update->ntimestep - (t0+t_equil); 

  if( (t >= 0) && (t <= t_switch) ) {
    lambda  =  switch_func(t/t_switch);
    dlambda = dswitch_func(t/t_switch);
  }

  if( (t >= t_equil+t_switch) && (t <= (t_equil+2*t_switch)) ) {
    lambda  =    switch_func(1.0 - (t - t_switch - t_equil)/t_switch);
    dlambda = - dswitch_func(1.0 - (t - t_switch - t_equil)/t_switch);
  }
}

/* ---------------------------------------------------------------------- */

double FixTIRS::compute_vector(int n)
{
  linfo[0] = lambda;
  linfo[1] = dlambda;
  return linfo[n];
}

/* ---------------------------------------------------------------------- */

double FixTIRS::switch_func(double t)
{
  if (sf == 2) return l_initial / (1 + t * (l_initial/l_final - 1));
  if (sf == 3) return l_initial / (1 + log2(1+t) * (l_initial/l_final - 1));
  
  // Default option is sf = 1.
  return l_initial + (l_final - l_initial) * t;
}

/* ---------------------------------------------------------------------- */

double FixTIRS::dswitch_func(double t)
{
  double aux = (1.0/l_initial - 1.0/l_final);
  if (sf == 2) return lambda * lambda * aux / t_switch;
  if (sf == 3) return lambda * lambda * aux / (t_switch * log(2) * (1 + t));

  // Default option is sf = 1.
  return (l_final-l_initial)/t_switch;
}
