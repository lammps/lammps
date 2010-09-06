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

#include "stdlib.h"
#include "string.h"
#include "fix_wall_reflect.h"
#include "atom.h"
#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixWallReflect::FixWallReflect(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix wall/reflect command");

  xloflag = xhiflag = yloflag = yhiflag = zloflag = zhiflag = NONE;
  xlostr = xhistr = ylostr = yhistr = zlostr = zhistr = NULL;
  int scaleflag = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"xlo") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	xloflag = EDGE;
	xlo = domain->boxlo[0];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	xloflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	xlostr = new char[n];
	strcpy(xlostr,&arg[iarg+1][2]);
      } else {
	xloflag = CONSTANT;
	xlo = atof(arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"xhi") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	xhiflag = EDGE;
	xhi = domain->boxhi[0];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	xhiflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	xhistr = new char[n];
	strcpy(xhistr,&arg[iarg+1][2]);
      } else {
	xhiflag = CONSTANT;
	xhi = atof(arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ylo") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	yloflag = EDGE;
	ylo = domain->boxlo[1];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	yloflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	ylostr = new char[n];
	strcpy(ylostr,&arg[iarg+1][2]);
      } else {
	yloflag = CONSTANT;
	ylo = atof(arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"yhi") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	yhiflag = EDGE;
	yhi = domain->boxhi[1];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	yhiflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	yhistr = new char[n];
	strcpy(yhistr,&arg[iarg+1][2]);
      } else {
	yhiflag = CONSTANT;
	yhi = atof(arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"zlo") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	zloflag = EDGE;
	zlo = domain->boxlo[2];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	zloflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	zlostr = new char[n];
	strcpy(zlostr,&arg[iarg+1][2]);
      } else {
	zloflag = CONSTANT;
	zlo = atof(arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"zhi") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix wall/reflect command");
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
	zhiflag = EDGE;
	zhi = domain->boxhi[2];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	zhiflag = VARIABLE;
	int n = strlen(&arg[iarg+1][2]) + 1;
	zhistr = new char[n];
	strcpy(zhistr,&arg[iarg+1][2]);
      } else {
	zhiflag = CONSTANT;
	zhi = atof(arg[iarg+1]);
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal wall/reflect command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal fix wall/reflect command");
      iarg += 2;
    } else error->all("Illegal fix wall/reflect command");
  }

  // setup scaling if any wall positions were set as CONSTANT
  // apply scaling factors to wall positions

  if (xloflag == CONSTANT || xhiflag == CONSTANT ||
      yloflag == CONSTANT || yhiflag == CONSTANT ||
      zloflag == CONSTANT || zhiflag == CONSTANT) {

    if (scaleflag && domain->lattice == NULL)
      error->all("Use of fix wall/reflect with undefined lattice");

    double xscale,yscale,zscale;
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    if (xloflag == CONSTANT) xlo *= xscale;
    if (xhiflag == CONSTANT) xhi *= xscale;
    if (yloflag == CONSTANT) ylo *= xscale;
    if (yhiflag == CONSTANT) yhi *= xscale;
    if (zloflag == CONSTANT) zlo *= xscale;
    if (zhiflag == CONSTANT) zhi *= xscale;
  }

  // error checks

  if (!xloflag && !xhiflag && !yloflag && !yhiflag && !zloflag && !zhiflag)
    error->all("Illegal fix wall/reflect command");

  if ((xloflag || xhiflag) && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if ((yloflag || yhiflag) && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if ((zloflag || zhiflag) && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
}

/* ---------------------------------------------------------------------- */

FixWallReflect::~FixWallReflect()
{
  delete [] xlostr;
  delete [] xhistr;
  delete [] ylostr;
  delete [] yhistr;
  delete [] zlostr;
  delete [] zhistr;
}

/* ---------------------------------------------------------------------- */

int FixWallReflect::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::init()
{
  if (xlostr) {
    xlovar = input->variable->find(xlostr);
    if (xlovar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(xlovar))
      error->all("Variable for fix wall/reflect is invalid style");
  }
  if (xhistr) {
    xhivar = input->variable->find(xhistr);
    if (xhivar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(xhivar))
      error->all("Variable for fix wall/reflect is invalid style");
  }
  if (ylostr) {
    ylovar = input->variable->find(ylostr);
    if (ylovar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(ylovar))
      error->all("Variable for fix wall/reflect is invalid style");
  }
  if (yhistr) {
    yhivar = input->variable->find(yhistr);
    if (yhivar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(yhivar))
      error->all("Variable for fix wall/reflect is invalid style");
  }
  if (zlostr) {
    zlovar = input->variable->find(zlostr);
    if (zlovar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(zlovar))
      error->all("Variable for fix wall/reflect is invalid style");
  }
  if (zhistr) {
    zhivar = input->variable->find(zhistr);
    if (zhivar < 0)
      error->all("Variable name for fix wall/relect does not exist");
    if (!input->variable->equalstyle(zhivar))
      error->all("Variable for fix wall/reflect is invalid style");
  }

  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0) 
    error->warning("Should not allow rigid bodies to bounce off "
		   "relecting walls");
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::post_integrate()
{
  if (xloflag == VARIABLE) xlo = input->variable->compute_equal(xlovar);
  if (xhiflag == VARIABLE) xhi = input->variable->compute_equal(xhivar);
  if (yloflag == VARIABLE) ylo = input->variable->compute_equal(ylovar);
  if (yhiflag == VARIABLE) yhi = input->variable->compute_equal(yhivar);
  if (zloflag == VARIABLE) zlo = input->variable->compute_equal(zlovar);
  if (zhiflag == VARIABLE) zhi = input->variable->compute_equal(zhivar);

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (xloflag && x[i][0] < xlo) {
	x[i][0] = xlo + (xlo - x[i][0]);
	v[i][0] = -v[i][0];
      }
      if (xhiflag && x[i][0] > xhi) {
	x[i][0] = xhi - (x[i][0] - xhi);
	v[i][0] = -v[i][0];
      }
      if (yloflag && x[i][1] < ylo) {
	x[i][1] = ylo + (ylo - x[i][1]);
	v[i][1] = -v[i][1];
      }
      if (yhiflag && x[i][1] > yhi) {
	x[i][1] = yhi - (x[i][1] - yhi);
	v[i][1] = -v[i][1];
      }
      if (zloflag && x[i][2] < zlo) {
	x[i][2] = zlo + (zlo - x[i][2]);
	v[i][2] = -v[i][2];
      }
      if (zhiflag && x[i][2] > zhi) {
	x[i][2] = zhi - (x[i][2] - zhi);
	v[i][2] = -v[i][2];
      }
    }
  }
}
