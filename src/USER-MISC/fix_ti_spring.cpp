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
#include "fix_ti_spring.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS; 
using namespace FixConst;  

/* ---------------------------------------------------------------------- */

FixTISpring::FixTISpring(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {

  if (narg < 6 || narg > 8)
    error->all(FLERR,"Illegal fix ti/spring command");

  // Flags.
  restart_peratom = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // Spring constant.
  k = atof(arg[3]);
  espring = 0.0;
  if (k <= 0.0) error->all(FLERR,"Illegal fix ti/spring command");

  // Initial position.
  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i], image[i], xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // Time variables.
  t_switch = atof(arg[4]); // Switching time.
  t_equil = atof(arg[5]);  // Equilibration time.
  t0 = update->ntimestep;  // Initial time.
  if (t_switch < 0.0) error->all(FLERR,"Illegal fix ti/spring command");
  if (t_equil  < 0.0) error->all(FLERR,"Illegal fix ti/spring command");

  // Coupling parameter initialization.
  sf = 1;
  if (narg > 6) {
    if (strcmp(arg[6], "function") == 0) sf = atoi(arg[7]);
    else error->all(FLERR,"Illegal fix ti/spring switching function");
    if ((sf!=1) && (sf!=2)) 
      error->all(FLERR,"Illegal fix ti/spring switching function");
  }
  lambda  =  switch_func(0); 
  dlambda = dswitch_func(0); 
}

/* ---------------------------------------------------------------------- */

FixTISpring::~FixTISpring()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored array
  memory->destroy(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixTISpring::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::setup(int vflag)
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

void FixTISpring::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::post_force(int vflag)
{

  // If on the first equilibration do not calculate forces.
  int t = update->ntimestep - t0;
  espring = 0.0;
  if(t < t_equil) return;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - xoriginal[i][0];
      dy = unwrap[1] - xoriginal[i][1];
      dz = unwrap[2] - xoriginal[i][2];
      f[i][0] = (1-lambda) * f[i][0] + lambda * (-k*dx);
      f[i][1] = (1-lambda) * f[i][1] + lambda * (-k*dy);
      f[i][2] = (1-lambda) * f[i][2] + lambda * (-k*dz);
      espring += k * (dx*dx + dy*dy + dz*dz);
    }

  espring *= 0.5;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::initial_integrate(int vflag) 
{ 
  // Update the coupling parameter value.
  double t = update->ntimestep - (t0+t_equil); 

  if( (t >= 0) && (t <= t_switch) ) {
    lambda  =  switch_func(t/t_switch); 
    dlambda = dswitch_func(t/t_switch); 
  }

  if( (t >= t_equil+t_switch) && (t <= (t_equil+2*t_switch)) ) {
    lambda  =    switch_func(1.0 - (t - t_switch - t_equil)/t_switch ); 
    dlambda = - dswitch_func(1.0 - (t - t_switch - t_equil)/t_switch ); 
  }
} 

/* ----------------------------------------------------------------------
   energy of stretched springs
------------------------------------------------------------------------- */

double FixTISpring::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   information about coupling parameter
------------------------------------------------------------------------- */

double FixTISpring::compute_vector(int n)
{
  linfo[0] = lambda;
  linfo[1] = dlambda;
  return linfo[n];
}

/* ----------------------------------------------------------------------
     memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTISpring::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
     allocate atom-based array
------------------------------------------------------------------------- */

void FixTISpring::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"fix_ti/spring:xoriginal");
}

/* ----------------------------------------------------------------------
     copy values within local atom-based array
------------------------------------------------------------------------- */

void FixTISpring::copy_arrays(int i, int j)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
    pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixTISpring::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
    unpack values in local atom-based array from exchange with another proc
 ------------------------------------------------------------------------- */

int FixTISpring::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
    pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTISpring::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}


/* ----------------------------------------------------------------------
    unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTISpring::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;
  
  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
     maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTISpring::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
     size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTISpring::size_restart(int nlocal)
{
  return 4;
}

/* ----------------------------------------------------------------------
     Switching function.
------------------------------------------------------------------------- */

double FixTISpring::switch_func(double t) 
{
  if (sf == 1) return t;

  double t2 = t*t;
  double t5 = t2*t2*t;
  return ((70.0*t2*t2 - 315.0*t2*t + 540.0*t2 - 420.0*t + 126.0)*t5);
}

/* ----------------------------------------------------------------------
     Switching function derivative.
------------------------------------------------------------------------- */

double FixTISpring::dswitch_func(double t) 
{
  if(sf == 1) return 1.0/t_switch;

  double t2 = t*t;
  double t4 = t2*t2;
  return ((630*t2*t2 - 2520*t2*t + 3780*t2 - 2520*t + 630)*t4) / t_switch;
}
