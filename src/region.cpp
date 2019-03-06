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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "region.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "math_extra.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Region::Region(LAMMPS *lmp, int /*narg*/, char **arg) :
  Pointers(lmp),
  id(NULL), style(NULL), contact(NULL), list(NULL),
  xstr(NULL), ystr(NULL), zstr(NULL), tstr(NULL)
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

  size_restart = 5;
  reset_vel();
  copymode = 0;
  list = NULL;
  nregion = 1;

  extent_xlo = 0.0;
  extent_ylo = 0.0;
  extent_zlo = 0.0;
  extent_xhi = 0.0;
  extent_yhi = 0.0;
  extent_zhi = 0.0;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  if (copymode) return;

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
  vel_timestep = -1;
}

/* ----------------------------------------------------------------------
   return 1 if region is dynamic (moves/rotates) or has variable shape
   else return 0 if static
------------------------------------------------------------------------- */

int Region::dynamic_check()
{
  if (dynamic || varshape) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   called before looping over atoms with match() or surface()
   this insures any variables used by region are invoked once per timestep
     also insures variables are invoked by all procs even those w/out atoms
     necessary if equal-style variable invokes global operation
   with MPI_Allreduce, e.g. xcm() or count()
------------------------------------------------------------------------- */

void Region::prematch()
{
  if (varshape) shape_update();
  if (dynamic) pretransform();
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
  if (dynamic) inverse_transform(x,y,z);
  if (openflag) return 1;
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

  if (dynamic) {
    xorig[0] = x;
    xorig[1] = y;
    xorig[2] = z;
    inverse_transform(x,y,z);
  }

  xnear[0] = x;
  xnear[1] = y;
  xnear[2] = z;

  if (!openflag) {
    if (interior) ncontact = surface_interior(xnear,cutoff);
    else ncontact = surface_exterior(xnear,cutoff);
  } else {
    // one of surface_int/ext() will return 0
    // so no need to worry about offset of contact indices
    ncontact = surface_exterior(xnear,cutoff) + surface_interior(xnear,cutoff);
  }

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
  contact[n].radius = 0;
  contact[n].delx = delx;
  contact[n].dely = dely;
  contact[n].delz = delz;
}

/* ----------------------------------------------------------------------
   pre-compute dx,dy,dz and theta for a moving/rotating region
   called once for the region before per-atom loop, via prematch()
------------------------------------------------------------------------- */

void Region::pretransform()
{
  if (moveflag) {
    if (xstr) dx = input->variable->compute_equal(xvar);
    if (ystr) dy = input->variable->compute_equal(yvar);
    if (zstr) dz = input->variable->compute_equal(zvar);
  }
  if (rotateflag) theta = input->variable->compute_equal(tvar);
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in region space to moved space
   rotate first (around original P), then displace
------------------------------------------------------------------------- */

void Region::forward_transform(double &x, double &y, double &z)
{
  if (rotateflag) rotate(x,y,z,theta);
  if (moveflag) {
    x += dx;
    y += dy;
    z += dz;
  }
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

void Region::inverse_transform(double &x, double &y, double &z)
{
  if (moveflag) {
    x -= dx;
    y -= dy;
    z -= dz;
  }
  if (rotateflag) rotate(x,y,z,-theta);
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   R = vector axis of rotation
   P = point = point to rotate around
   R0 = runit = unit vector for R
   X0 = x,y,z = initial coord of atom
   D = X0 - P = vector from P to X0
   C = (D dot R0) R0 = projection of D onto R, i.e. Dparallel
   A = D - C = vector from R line to X0, i.e. Dperp
   B = R0 cross A = vector perp to A in plane of rotation, same len as A
   A,B define plane of circular rotation around R line
   new x,y,z = P + C + A cos(angle) + B sin(angle)
------------------------------------------------------------------------- */

void Region::rotate(double &x, double &y, double &z, double angle)
{
  double a[3],b[3],c[3],d[3],disp[3];

  double sine = sin(angle);
  double cosine = cos(angle);
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  double x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
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

  openflag = 0;
  for (int i  = 0; i < 6; i++) open_faces[i] = 0;

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

    } else if (strcmp(arg[iarg],"open") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal region command");
      int iface = force->inumeric(FLERR,arg[iarg+1]);
      if (iface < 1 || iface > 6) error->all(FLERR,"Illegal region command");
      // additional checks on valid face index are done by region classes
      open_faces[iface-1] = 1;
      openflag = 1;
      iarg += 2;
    }
    else error->all(FLERR,"Illegal region command");
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

/* ----------------------------------------------------------------------
   find nearest point to C on line segment A,B and return it as D
   project (C-A) onto (B-A)
   t = length of that projection, normalized by length of (B-A)
   t <= 0, C is closest to A
   t >= 1, C is closest to B
   else closest point is between A and B
------------------------------------------------------------------------- */

void Region::point_on_line_segment(double *a, double *b,
                                   double *c, double *d)
{
  double ba[3],ca[3];

  MathExtra::sub3(b,a,ba);
  MathExtra::sub3(c,a,ca);
  double t = MathExtra::dot3(ca,ba) / MathExtra::dot3(ba,ba);
  if (t <= 0.0) {
    d[0] = a[0];
    d[1] = a[1];
    d[2] = a[2];
  } else if (t >= 1.0) {
    d[0] = b[0];
    d[1] = b[1];
    d[2] = b[2];
  } else {
    d[0] = a[0] + t*ba[0];
    d[1] = a[1] + t*ba[1];
    d[2] = a[2] + t*ba[2];
  }
}

/* ----------------------------------------------------------------------
   infer translational and angular velocity of region
   necessary b/c motion variables are for displacement & theta
     there is no analytic formula for v & omega
   prev[4] contains values of dx,dy,dz,theta at previous step
     used for difference, then updated to current step values
   dt is time elapsed since previous step
   rpoint = point updated by current displacement
   called by fix wall/gran/region every timestep
------------------------------------------------------------------------- */

void Region::set_velocity()
{
  if (vel_timestep == update->ntimestep) return;
  vel_timestep = update->ntimestep;
  if (moveflag) {
    if (update->ntimestep > 0) {
      v[0] = (dx - prev[0])/update->dt;
      v[1] = (dy - prev[1])/update->dt;
      v[2] = (dz - prev[2])/update->dt;
    }
    else v[0] = v[1] = v[2] = 0.0;
    prev[0] = dx;
    prev[1] = dy;
    prev[2] = dz;
  }

  if (rotateflag) {
    rpoint[0] = point[0] + dx;
    rpoint[1] = point[1] + dy;
    rpoint[2] = point[2] + dz;
    if (update->ntimestep > 0) {
      double angvel = (theta-prev[3]) / update->dt;
      omega[0] = angvel*axis[0];
      omega[1] = angvel*axis[1];
      omega[2] = angvel*axis[2];
    }
    else omega[0] = omega[1] = omega[2] = 0.0;
    prev[3] = theta;
  }

  if (varshape){
    set_velocity_shape();
  }
}

/* ----------------------------------------------------------------------
   compute velocity of wall for given contact
   since contacts only store delx/y/z, need to pass particle coords
     to compute contact point
   called by fix/wall/gran/region every contact every timestep
------------------------------------------------------------------------- */

void Region::velocity_contact(double *vwall, double *x, int ic)
{
  double xc[3];

  vwall[0] = vwall[1] = vwall[2] = 0.0;

  if (moveflag){
    vwall[0] = v[0];
    vwall[1] = v[1];
    vwall[2] = v[2];
  }
  if (rotateflag){
    xc[0] = x[0] - contact[ic].delx;
    xc[1] = x[1] - contact[ic].dely;
    xc[2] = x[2] - contact[ic].delz;
    vwall[0] += omega[1]*(xc[2] - rpoint[2]) - omega[2]*(xc[1] - rpoint[1]);
    vwall[1] += omega[2]*(xc[0] - rpoint[0]) - omega[0]*(xc[2] - rpoint[2]);
    vwall[2] += omega[0]*(xc[1] - rpoint[1]) - omega[1]*(xc[0] - rpoint[0]);
  }

  if (varshape && contact[ic].varflag) velocity_contact_shape(vwall, xc);
}


/* ----------------------------------------------------------------------
   increment length of restart buffer based on region info
   used by restart of fix/wall/gran/region
------------------------------------------------------------------------- */

void Region::length_restart_string(int &n)
{
  n += sizeof(int) + strlen(id)+1 +
    sizeof(int) + strlen(style)+1 + sizeof(int) +
    size_restart*sizeof(double);
}

/* ----------------------------------------------------------------------
   region writes its current style, id, number of sub-regions, position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

void Region::write_restart(FILE *fp)
{
  int sizeid = (strlen(id)+1);
  int sizestyle = (strlen(style)+1);
  fwrite(&sizeid, sizeof(int), 1, fp);
  fwrite(id,1,sizeid,fp);
  fwrite(&sizestyle,sizeof(int),1,fp);
  fwrite(style,1,sizestyle,fp);
  fwrite(&nregion,sizeof(int),1,fp);
  fwrite(prev,sizeof(double),size_restart,fp);
}

/* ----------------------------------------------------------------------
   region reads style, id, number of sub-regions from restart file
   if they match current region, also read previous position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

int Region::restart(char *buf, int &n)
{
  int size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n],id) != 0)) return 0;
  n += size;

  size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n],style) != 0)) return 0;
  n += size;

  int restart_nreg = *((int *) (&buf[n]));
  n += sizeof(int);
  if (restart_nreg != nregion) return 0;

  memcpy(prev,&buf[n],size_restart*sizeof(double));
  return 1;
}

/* ----------------------------------------------------------------------
   set prev vector to zero
------------------------------------------------------------------------- */

void Region::reset_vel()
{
  for (int i = 0; i < size_restart; i++) prev[i] = 0;
}

void Region::get_bounding_box(double * xlo, double * xhi) const
{
  if(bboxflag) {
    xlo[0] = extent_xlo;
    xlo[1] = extent_ylo;
    xlo[2] = extent_zlo;
    xhi[0] = extent_xhi;
    xhi[1] = extent_yhi;
    xhi[2] = extent_zhi;
  } else {
    xlo[0] = domain->boxlo[0];
    xlo[1] = domain->boxlo[1];
    xlo[2] = domain->boxlo[2];
    xhi[0] = domain->boxhi[0];
    xhi[1] = domain->boxhi[1];
    xhi[2] = domain->boxhi[2];
  }
}
