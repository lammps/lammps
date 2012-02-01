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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{LJ93,LJ126,COLLOID,HARMONIC};

/* ---------------------------------------------------------------------- */

FixWallRegion::FixWallRegion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix wall/region command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // parse args

  iregion = domain->find_region(arg[3]);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/region does not exist");
  int n = strlen(arg[3]) + 1;
  idregion = new char[n];
  strcpy(idregion,arg[3]);

  if (strcmp(arg[4],"lj93") == 0) style = LJ93;
  else if (strcmp(arg[4],"lj126") == 0) style = LJ126;
  else if (strcmp(arg[4],"colloid") == 0) style = COLLOID;
  else if (strcmp(arg[4],"harmonic") == 0) style = HARMONIC;
  else error->all(FLERR,"Illegal fix wall/region command");

  epsilon = atof(arg[5]);
  sigma = atof(arg[6]);
  cutoff = atof(arg[7]);

  if (cutoff <= 0.0) error->all(FLERR,"Fix wall/region cutoff <= 0.0");

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWallRegion::~FixWallRegion()
{
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixWallRegion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallRegion::init()
{
  // set index and check validity of region

  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/region does not exist");

  // error checks for style COLLOID
  // insure all particles in group are extended particles

  if (style == COLLOID) {
    if (!atom->sphere_flag) 
      error->all(FLERR,"Fix wall/region colloid requires atom style sphere");

    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	if (radius[i] == 0.0) flag = 1;
    
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) 
      error->all(FLERR,"Fix wall/region colloid requires extended particles");
  }

  // setup coefficients for each style

  if (style == LJ93) {
    coeff1 = 6.0/5.0 * epsilon * pow(sigma,9.0);
    coeff2 = 3.0 * epsilon * pow(sigma,3.0);
    coeff3 = 2.0/15.0 * epsilon * pow(sigma,9.0);
    coeff4 = epsilon * pow(sigma,3.0);
    double rinv = 1.0/cutoff;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;
  } else if (style == LJ126) {
    coeff1 = 48.0 * epsilon * pow(sigma,12.0);
    coeff2 = 24.0 * epsilon * pow(sigma,6.0);
    coeff3 = 4.0 * epsilon * pow(sigma,12.0);
    coeff4 = 4.0 * epsilon * pow(sigma,6.0);
    double r2inv = 1.0/(cutoff*cutoff);
    double r6inv = r2inv*r2inv*r2inv;
    offset = r6inv*(coeff3*r6inv - coeff4);
  } else if (style == COLLOID) {
    coeff1 = -4.0/315.0 * epsilon * pow(sigma,6.0);
    coeff2 = -2.0/3.0 * epsilon;
    coeff3 = epsilon * pow(sigma,6.0)/7560.0;
    coeff4 = epsilon/6.0;
    double rinv = 1.0/cutoff;
    double r2inv = rinv*rinv;
    double r4inv = r2inv*r2inv;
    offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallRegion::setup(int vflag)
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

void FixWallRegion::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegion::post_force(int vflag)
{
  int i,m,n;
  double rinv,fx,fy,fz,tooclose;

  eflag = 0;
  ewall[0] = ewall[1] = ewall[2] = ewall[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  Region *region = domain->regions[iregion];
  int onflag = 0;

  // region->match() insures particle is in region or on surface, else error
  // if returned contact dist r = 0, is on surface, also an error
  // in COLLOID case, r <= radius is an error
  
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (!region->match(x[i][0],x[i][1],x[i][2])) {
	onflag = 1;
	continue;
      }
      if (style == COLLOID) tooclose = radius[i];
      else tooclose = 0.0;

      n = region->surface(x[i][0],x[i][1],x[i][2],cutoff);

      for (m = 0; m < n; m++) {
	if (region->contact[m].r <= tooclose) {
	  onflag = 1;
	  continue;
	} else rinv = 1.0/region->contact[m].r;

	if (style == LJ93) lj93(region->contact[m].r);
	else if (style == LJ126) lj126(region->contact[m].r);
	else if (style == COLLOID) colloid(region->contact[m].r,radius[i]);
	else harmonic(region->contact[m].r);

	ewall[0] += eng;
	fx = fwall * region->contact[m].delx * rinv;
	fy = fwall * region->contact[m].dely * rinv;
	fz = fwall * region->contact[m].delz * rinv;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;
	ewall[1] -= fx;
	ewall[2] -= fy;
	ewall[3] -= fz;
      }
    }

  if (onflag) error->one(FLERR,"Particle on or inside surface of region "
			 "used in fix wall/region");
}

/* ---------------------------------------------------------------------- */

void FixWallRegion::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallRegion::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallRegion::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallRegion::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,4,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}

/* ----------------------------------------------------------------------
   LJ 9/3 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

void FixWallRegion::lj93(double r)
{
  double rinv = 1.0/r;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  double r10inv = r4inv*r4inv*r2inv;
  fwall = coeff1*r10inv - coeff2*r4inv;
  eng = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv - offset;
}

/* ----------------------------------------------------------------------
   LJ 12/6 interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

void FixWallRegion::lj126(double r)
{
  double rinv = 1.0/r;
  double r2inv = rinv*rinv;
  double r6inv = r2inv*r2inv*r2inv;
  fwall = r6inv*(coeff1*r6inv - coeff2) * rinv;
  eng = r6inv*(coeff3*r6inv - coeff4) - offset;
}

/* ----------------------------------------------------------------------
   colloid interaction for finite-size particle of rad with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

void FixWallRegion::colloid(double r, double rad)
{
  double new_coeff2 = coeff2*rad*rad*rad;
  double diam = 2.0*rad;

  double rad2 = rad*rad;
  double rad4 = rad2*rad2;
  double rad8 = rad4*rad4;
  double delta2 = rad2 - r*r;
  double rinv = 1.0/delta2;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  double r8inv = r4inv*r4inv;
  fwall = coeff1*(rad8*rad + 27.0*rad4*rad2*rad*pow(r,2.0)
		  + 63.0*rad4*rad*pow(r,4.0)
		  + 21.0*rad2*rad*pow(r,6.0))*r8inv - new_coeff2*r2inv;

  double r2 = 0.5*diam - r;
  double rinv2 = 1.0/r2;
  double r2inv2 = rinv2*rinv2;
  double r4inv2 = r2inv2*r2inv2;
  double r3 = r + 0.5*diam;
  double rinv3 = 1.0/r3;
  double r2inv3 = rinv3*rinv3;
  double r4inv3 = r2inv3*r2inv3;
  eng = coeff3*((-3.5*diam+r)*r4inv2*r2inv2*rinv2
		+ (3.5*diam+r)*r4inv3*r2inv3*rinv3) -
    coeff4*((-diam*r+r2*r3*(log(-r2)-log(r3)))*
	    (-rinv2)*rinv3) - offset;
}

/* ----------------------------------------------------------------------
   harmonic interaction for particle with wall
   compute eng and fwall = magnitude of wall force
------------------------------------------------------------------------- */

void FixWallRegion::harmonic(double r)
{
  double dr = cutoff - r;
  fwall = 2.0*epsilon*dr;
  eng = epsilon*dr*dr;
}
