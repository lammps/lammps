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
#include "region.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Region::Region(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  varshape = 0;
  xstr = ystr = zstr = tstr = NULL;
  dx = dy = dz = 0.0;
  lastshape = lastdynamic = -1;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  delete [] id;
  delete [] style;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] tstr;
}

/* ---------------------------------------------------------------------- */

void Region::init()
{
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for region is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR,"Variable name for region does not exist");
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR,"Variable for region is not equal style");
  }
}

/* ----------------------------------------------------------------------
   return 1 if region is dynamic, 0 if static
   only primitive regions define it here
   union/intersect regions have their own dynamic_check()
------------------------------------------------------------------------- */

int Region::dynamic_check()
{
  return dynamic;
}

/* ----------------------------------------------------------------------
   determine if point x,y,z is a match to region volume
   XOR computes 0 if 2 args are the same, 1 if different
   note that inside() returns 1 for points on surface of region
   thus point on surface of exterior region will not match
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic, apply inverse transform to x,y,z
     unmove first, then unrotate, so don't have to change rotation point
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

int Region::match(double x, double y, double z)
{
  if (varshape && update->ntimestep != lastshape) {
    shape_update();
    lastshape = update->ntimestep;
  }

  if (dynamic) inverse_transform(x,y,z);

  return !(inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   generate list of contact points for interior or exterior regions
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic:
     before: inverse transform x,y,z (unmove, then unrotate)
     after: forward transform contact point xs,yx,zs (rotate, then move),
            then reset contact delx,dely,delz based on new contact point
            no need to do this if no rotation since delxyz doesn't change
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

int Region::surface(double x, double y, double z, double cutoff)
{
  int ncontact;
  double xs,ys,zs;
  double xnear[3],xorig[3];

  if (varshape && update->ntimestep != lastshape) {
    shape_update();
    lastshape = update->ntimestep;
  }

  if (dynamic) {
    xorig[0] = x;
    xorig[1] = y;
    xorig[2] = z;
    inverse_transform(x,y,z);
  }

  xnear[0] = x;
  xnear[1] = y;
  xnear[2] = z;

  if (interior) ncontact = surface_interior(xnear,cutoff);
  else ncontact = surface_exterior(xnear,cutoff);

  if (rotateflag && ncontact) {
    for (int i = 0; i < ncontact; i++) {
      xs = xnear[0] - contact[i].delx;
      ys = xnear[1] - contact[i].dely;
      zs = xnear[2] - contact[i].delz;
      forward_transform(xs,ys,zs);
      contact[i].delx = xorig[0] - xs;
      contact[i].dely = xorig[1] - ys;
      contact[i].delz = xorig[2] - zs;
    }
  }

  return ncontact;
}

/* ----------------------------------------------------------------------
   add a single contact at Nth location in contact array
   x = particle position
   xp,yp,zp = region surface point
------------------------------------------------------------------------- */

void Region::add_contact(int n, double *x, double xp, double yp, double zp)
{
  double delx = x[0] - xp;
  double dely = x[1] - yp;
  double delz = x[2] - zp;
  contact[n].r = sqrt(delx*delx + dely*dely + delz*delz);
  contact[n].delx = delx;
  contact[n].dely = dely;
  contact[n].delz = delz;
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in region space to moved space
   rotate first (around original P), then displace
------------------------------------------------------------------------- */

void Region::forward_transform(double &x, double &y, double &z)
{
  if (rotateflag) {
    if (update->ntimestep != lastdynamic)
      theta = input->variable->compute_equal(tvar);
    rotate(x,y,z,theta);
  }

  if (moveflag) {
    if (update->ntimestep != lastdynamic) {
      if (xstr) dx = input->variable->compute_equal(xvar);
      if (ystr) dy = input->variable->compute_equal(yvar);
      if (zstr) dz = input->variable->compute_equal(zvar);
    }
    x += dx;
    y += dy;
    z += dz;
  }

  lastdynamic = update->ntimestep;
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

void Region::inverse_transform(double &x, double &y, double &z)
{
  if (moveflag) {
    if (update->ntimestep != lastdynamic) {
      if (xstr) dx = input->variable->compute_equal(xvar);
      if (ystr) dy = input->variable->compute_equal(yvar);
      if (zstr) dz = input->variable->compute_equal(zvar);
    }
    x -= dx;
    y -= dy;
    z -= dz;
  }

  if (rotateflag) {
    if (update->ntimestep != lastdynamic)
      theta = input->variable->compute_equal(tvar);
    rotate(x,y,z,-theta);
  }

  lastdynamic = update->ntimestep;
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   P = point = vector = point of rotation
   R = vector = axis of rotation
   w = omega of rotation (from period)
   X0 = x,y,z = initial coord of atom
   R0 = runit = unit vector for R
   C = (X0 dot R0) R0 = projection of atom coord onto R
   D = X0 - P = vector from P to X0
   A = D - C = vector from R line to X0
   B = R0 cross A = vector perp to A in plane of rotation
   A,B define plane of circular rotation around R line
   x,y,z = P + C + A cos(w*dt) + B sin(w*dt)
------------------------------------------------------------------------- */

void Region::rotate(double &x, double &y, double &z, double angle)
{
  double a[3],b[3],c[3],d[3],disp[3];

  double sine = sin(angle);
  double cosine = cos(angle);
  double x0dotr = x*runit[0] + y*runit[1] + z*runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1]*a[2] - runit[2]*a[1];
  b[1] = runit[2]*a[0] - runit[0]*a[2];
  b[2] = runit[0]*a[1] - runit[1]*a[0];
  disp[0] = a[0]*cosine  + b[0]*sine;
  disp[1] = a[1]*cosine  + b[1]*sine;
  disp[2] = a[2]*cosine  + b[2]*sine;
  x = point[0] + c[0] + disp[0];
  y = point[1] + c[1] + disp[1];
  z = point[2] + c[2] + disp[2];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal region command");

  // option defaults

  interior = 1;
  scaleflag = 1;
  moveflag = rotateflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all(FLERR,"Illegal region command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"move") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal region command");
      if (strcmp(arg[iarg+1],"NULL") != 0) {
        if (strstr(arg[iarg+1],"v_") != arg[iarg+1])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr = new char[n];
        strcpy(xstr,&arg[iarg+1][2]);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0) {
        if (strstr(arg[iarg+2],"v_") != arg[iarg+2])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+2][2]) + 1;
        ystr = new char[n];
        strcpy(ystr,&arg[iarg+2][2]);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0) {
        if (strstr(arg[iarg+3],"v_") != arg[iarg+3])
          error->all(FLERR,"Illegal region command");
        int n = strlen(&arg[iarg+3][2]) + 1;
        zstr = new char[n];
        strcpy(zstr,&arg[iarg+3][2]);
      }
      moveflag = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+8 > narg) error->all(FLERR,"Illegal region command");
      if (strstr(arg[iarg+1],"v_") != arg[iarg+1])
        error->all(FLERR,"Illegal region command");
      int n = strlen(&arg[iarg+1][2]) + 1;
      tstr = new char[n];
      strcpy(tstr,&arg[iarg+1][2]);
      point[0] = force->numeric(FLERR,arg[iarg+2]);
      point[1] = force->numeric(FLERR,arg[iarg+3]);
      point[2] = force->numeric(FLERR,arg[iarg+4]);
      axis[0] = force->numeric(FLERR,arg[iarg+5]);
      axis[1] = force->numeric(FLERR,arg[iarg+6]);
      axis[2] = force->numeric(FLERR,arg[iarg+7]);
      rotateflag = 1;
      iarg += 8;
    } else error->all(FLERR,"Illegal region command");
  }

  // error check

  if ((moveflag || rotateflag) &&
      (strcmp(style,"union") == 0 || strcmp(style,"intersect") == 0))
    error->all(FLERR,"Region union or intersect cannot be dynamic");

  // setup scaling

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  if (rotateflag) {
    point[0] *= xscale;
    point[1] *= yscale;
    point[2] *= zscale;
  }

  // runit = unit vector along rotation axis

  if (rotateflag) {
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Region cannot have 0 length rotation vector");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }

  if (moveflag || rotateflag) dynamic = 1;
  else dynamic = 0;
}
