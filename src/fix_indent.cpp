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
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_indent.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,SPHERE,CYLINDER,PLANE};
enum{INSIDE,OUTSIDE};

/* ---------------------------------------------------------------------- */

FixIndent::FixIndent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix indent command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  scalar_vector_freq = 1;
  extscalar = 1;
  extvector = 1;

  k = atof(arg[3]);

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of fix indent with undefined lattice");

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors to force constant, velocity, and geometry

  k /= xscale;
  k3 = k/3.0;
  vx *= xscale;
  vy *= yscale;
  vz *= zscale;

  if (istyle == SPHERE) {
    x0 *= xscale;
    y0 *= yscale;
    z0 *= zscale;
    r0_stop *= xscale;
    r0_start *= xscale;
  } else if (istyle == CYLINDER) {
    if (cdim == 0) {
      c1 *= yscale;
      c2 *= zscale;
      r0_stop *= xscale;
      r0_start *= xscale;
    } else if (cdim == 1) {
      c1 *= xscale;
      c2 *= zscale;
      r0_stop *= yscale;
      r0_start *= yscale;
    } else if (cdim == 2) {
      c1 *= xscale;
      c2 *= yscale;
      r0_stop *= zscale;
      r0_start *= zscale;
    }
  } else if (istyle == PLANE) {
    if (cdim == 0) planepos *= xscale;
    else if (cdim == 1) planepos *= yscale;
    else if (cdim == 2) planepos *= zscale;
  } else error->all("Illegal fix indent command");

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixIndent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndent::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixIndent::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force(int vflag)
{
  // set current r0
  // for minimization, always set to r0_stop

  double r0;
  if (!radflag || update->whichflag == 2) r0 = r0_stop;
  else {
    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;
    r0 = r0_start + delta * (r0_stop-r0_start);
  }

  // indenter values, 0 = energy, 1-3 = force components

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // spherical indenter

  if (istyle == SPHERE) {

    // ctr = current indenter center from original x0,y0,z0
    // remap into periodic box

    double ctr[3];
    double delta = (update->ntimestep - update->beginstep) * update->dt;
    ctr[0] = x0 + delta*vx;
    ctr[1] = y0 + delta*vy;
    ctr[2] = z0 + delta*vz;
    domain->remap(ctr);

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag,fx,fy,fz;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	delx = x[i][0] - ctr[0];
	dely = x[i][1] - ctr[1];
	delz = x[i][2] - ctr[2];
	domain->minimum_image(delx,dely,delz);
	r = sqrt(delx*delx + dely*dely + delz*delz);
	if (side == OUTSIDE) {
	  dr = r - r0;
	  fmag = k*dr*dr;
	} else {
	  dr = r0 - r;
	  fmag = -k*dr*dr;
	}
	if (dr >= 0.0) continue;
	fx = delx*fmag/r;
	fy = dely*fmag/r;
	fz = delz*fmag/r;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;
	indenter[0] -= k3 * dr*dr*dr;
	indenter[1] -= fx;
	indenter[2] -= fy;
	indenter[3] -= fz;
      }

  // cylindrical indenter

  } else if (istyle == CYLINDER) {

    // ctr = current indenter axis from original c1,c2
    // remap into periodic box
    // 3rd coord is just near box for remap(), since isn't used
	      
    double ctr[3];
    double delta = (update->ntimestep - update->beginstep) * update->dt;
    if (cdim == 0) {
      ctr[0] = domain->boxlo[0];
      ctr[1] = c1 + delta*vy;
      ctr[2] = c2 + delta*vz;
    } else if (cdim == 1) {
      ctr[0] = c1 + delta*vx;
      ctr[1] = domain->boxlo[1];
      ctr[2] = c2 + delta*vz;
    } else {
      ctr[0] = c1 + delta*vx;
      ctr[1] = c2 + delta*vy;
      ctr[2] = domain->boxlo[2];
    }
    domain->remap(ctr);

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    double delx,dely,delz,r,dr,fmag,fx,fy,fz;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (cdim == 0) {
	  delx = 0;
	  dely = x[i][1] - ctr[1];
	  delz = x[i][2] - ctr[2];
	} else if (cdim == 1) {
	  delx = x[i][0] - ctr[0];
	  dely = 0;
	  delz = x[i][2] - ctr[2];
	} else {
	  delx = x[i][0] - ctr[0];
	  dely = x[i][1] - ctr[1];
	  delz = 0;
	}
	domain->minimum_image(delx,dely,delz);
	r = sqrt(delx*delx + dely*dely + delz*delz);
	if (side == OUTSIDE) {
	  dr = r - r0;
	  fmag = k*dr*dr;
	} else {
	  dr = r0 - r;
	  fmag = -k*dr*dr;
	}
	if (dr >= 0.0) continue;
	fx = delx*fmag/r;
	fy = dely*fmag/r;
	fz = delz*fmag/r;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;
	indenter[0] -= k3 * dr*dr*dr;
	indenter[1] -= fx;
	indenter[2] -= fy;
	indenter[3] -= fz;
      }

  // planar indenter

  } else {

    // posnew = current coord of plane from original planepos
	      
    double delta = (update->ntimestep - update->beginstep) * update->dt;
    double posnew;
    if (cdim == 0) posnew = planepos + delta*vx;
    else if (cdim == 1) posnew = planepos + delta*vy;
    else if (cdim == 2) posnew = planepos + delta*vz;
    
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    
    double dr,fatom;
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	dr = planeside * (posnew - x[i][cdim]);
	if (dr >= 0.0) continue;
	fatom = -planeside * k*dr*dr;
	f[i][cdim] += fatom;
	indenter[0] -= k3 * dr*dr*dr;
	indenter[cdim+1] -= fatom;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixIndent::compute_scalar()
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[0];
}

/* ----------------------------------------------------------------------
   components of force on indenter
------------------------------------------------------------------------- */

double FixIndent::compute_vector(int n)
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[n+1];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line 
------------------------------------------------------------------------- */

void FixIndent::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal fix indent command");

  istyle = NONE;
  vx = vy = vz = 0.0;
  radflag = 0;
  r0_start = 0.0;
  scaleflag = 1;
  side = OUTSIDE;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      x0 = atof(arg[iarg+1]);
      y0 = atof(arg[iarg+2]);
      z0 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = SPHERE;
      iarg += 5;
    } else if (strcmp(arg[iarg],"cylinder") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"x") == 0) cdim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) cdim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) cdim = 2;
      else error->all("Illegal fix indent command");
      c1 = atof(arg[iarg+2]);
      c2 = atof(arg[iarg+3]);
      r0_stop = atof(arg[iarg+4]);
      istyle = CYLINDER;
      iarg += 5;
    } else if (strcmp(arg[iarg],"plane") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"x") == 0) cdim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) cdim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) cdim = 2;
      else error->all("Illegal fix indent command");
      planepos = atof(arg[iarg+2]);
      if (strcmp(arg[iarg+3],"lo") == 0) planeside = -1;
      else if (strcmp(arg[iarg+3],"hi") == 0) planeside = 1;
      else error->all("Illegal fix indent command");
      istyle = PLANE;
      iarg += 4;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix indent command");
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rstart") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      radflag = 1;
      r0_start = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal fix indent command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix indent command");
      if (strcmp(arg[iarg+1],"in") == 0) side = INSIDE;
      else if (strcmp(arg[iarg+1],"out") == 0) side = OUTSIDE;
      else error->all("Illegal fix indent command");
      iarg += 2;
    } else error->all("Illegal fix indent command");
  }
}
