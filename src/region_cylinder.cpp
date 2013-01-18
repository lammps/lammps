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
#include "region_cylinder.h"
#include "update.h"
#include "domain.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20
enum{CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

RegCylinder::RegCylinder(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"x") && strcmp(arg[2],"y") && strcmp(arg[2],"z"))
    error->all(FLERR,"Illegal region cylinder command");
  axis = arg[2][0];

  if (axis == 'x') {
    c1 = yscale*atof(arg[3]);
    c2 = zscale*atof(arg[4]);
  } else if (axis == 'y') {
    c1 = xscale*atof(arg[3]);
    c2 = zscale*atof(arg[4]);
  } else if (axis == 'z') {
    c1 = xscale*atof(arg[3]);
    c2 = yscale*atof(arg[4]);
  }

  rstr = NULL;
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    rstr = new char[n];
    strcpy(rstr,&arg[5][2]);
    radius = 0.0;
    rstyle = VARIABLE;
    varshape = 1;
    variable_check();
    shape_update();
  } else {
    radius = atof(arg[5]);
    if (axis == 'x') radius *= xscale;
    else radius *= xscale;
    rstyle = CONSTANT;
  }

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[0];
      else lo = domain->boxlo_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[1];
      else lo = domain->boxlo_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[6],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[2];
      else lo = domain->boxlo_bound[2];
    }
  } else {
    if (axis == 'x') lo = xscale*atof(arg[6]);
    if (axis == 'y') lo = yscale*atof(arg[6]);
    if (axis == 'z') lo = zscale*atof(arg[6]);
  }

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[0];
      else hi = domain->boxhi_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[1];
      else hi = domain->boxhi_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[7],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[2];
      else hi = domain->boxhi_bound[2];
    }
  } else {
    if (axis == 'x') hi = xscale*atof(arg[7]);
    if (axis == 'y') hi = yscale*atof(arg[7]);
    if (axis == 'z') hi = zscale*atof(arg[7]);
  }

  // error check

  if (radius <= 0.0) error->all(FLERR,"Illegal region cylinder command");

  // extent of cylinder
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 - radius;
      extent_yhi = c1 + radius;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'y') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'z') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = c2 - radius;
      extent_yhi = c2 + radius;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else bboxflag = 0;

  // particle could be contact cylinder surface and 2 ends

  cmax = 3;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegCylinder::~RegCylinder()
{
  delete [] rstr;
  delete [] contact;
}

/* ---------------------------------------------------------------------- */

void RegCylinder::init()
{
  Region::init();
  if (rstr) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCylinder::inside(double x, double y, double z)
{
  double del1,del2,dist;
  int inside;

  if (axis == 'x') {
    del1 = y - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && x >= lo && x <= hi) inside = 1;
    else inside = 0;
  } else if (axis == 'y') {
    del1 = x - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && y >= lo && y <= hi) inside = 1;
    else inside = 0;
  } else {
    del1 = x - c1;
    del2 = y - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && z >= lo && z <= hi) inside = 1;
    else inside = 0;
  }

  return inside;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of cylinder
   can be one contact for each of 3 cylinder surfaces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
   special case: no contact with curved surf if x is on center axis
------------------------------------------------------------------------- */

int RegCylinder::surface_interior(double *x, double cutoff)
{
  double del1,del2,r,delta;

  int n = 0;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // x is exterior to cylinder

    if (r > radius || x[0] < lo || x[0] > hi) return 0;

    // x is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0) {
      contact[n].r = delta;
      contact[n].delx = 0.0;
      contact[n].dely = del1*(1.0-radius/r);
      contact[n].delz = del2*(1.0-radius/r);
      n++;
    }
    delta = x[0] - lo;
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delx = delta;
      contact[n].dely = contact[n].delz = 0.0;
      n++;
    }
    delta = hi - x[0];
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delx = -delta;
      contact[n].dely = contact[n].delz = 0.0;
      n++;
    }

  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // y is exterior to cylinder

    if (r > radius || x[1] < lo || x[1] > hi) return 0;

    // y is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0) {
      contact[n].r = delta;
      contact[n].delx = del1*(1.0-radius/r);
      contact[n].dely = 0.0;
      contact[n].delz = del2*(1.0-radius/r);
      n++;
    }
    delta = x[1] - lo;
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].dely = delta;
      contact[n].delx = contact[n].delz = 0.0;
      n++;
    }
    delta = hi - x[1];
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].dely = -delta;
      contact[n].delx = contact[n].delz = 0.0;
      n++;
    }

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // z is exterior to cylinder

    if (r > radius || x[2] < lo || x[2] > hi) return 0;

    // z is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0) {
      contact[n].r = delta;
      contact[n].delx = del1*(1.0-radius/r);
      contact[n].dely = del2*(1.0-radius/r);
      contact[n].delz = 0.0;
      n++;
    }
    delta = x[2] - lo;
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delz = delta;
      contact[n].delx = contact[n].dely = 0.0;
      n++;
    }
    delta = hi - x[2];
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delz = -delta;
      contact[n].delx = contact[n].dely = 0.0;
      n++;
    }
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of cylinder
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
------------------------------------------------------------------------- */

int RegCylinder::surface_exterior(double *x, double cutoff)
{
  double del1,del2,r;
  double xp,yp,zp;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // x is far enough from cylinder that there is no contact
    // x is interior to cylinder

    if (r >= radius+cutoff || x[0] <= lo-cutoff || x[0] >= hi+cutoff) return 0;
    if (r < radius && x[0] > lo && x[0] < hi) return 0;

    // x is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    if (r > radius) {
      yp = c1 + del1*radius/r;
      zp = c2 + del2*radius/r;
    } else {
      yp = x[1];
      zp = x[2];
    }
    if (x[0] < lo) xp = lo;
    else if (x[0] > hi) xp = hi;
    else xp = x[0];

    add_contact(0,x,xp,yp,zp);
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // y is far enough from cylinder that there is no contact
    // y is interior to cylinder

    if (r >= radius+cutoff || x[1] <= lo-cutoff || x[1] >= hi+cutoff) return 0;
    if (r < radius && x[1] > lo && x[1] < hi) return 0;

    // y is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    if (r > radius) {
      xp = c1 + del1*radius/r;
      zp = c2 + del2*radius/r;
    } else {
      xp = x[0];
      zp = x[2];
    }
    if (x[1] < lo) yp = lo;
    else if (x[1] > hi) yp = hi;
    else yp = x[1];

    add_contact(0,x,xp,yp,zp);
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1*del1 + del2*del2);

    // z is far enough from cylinder that there is no contact
    // z is interior to cylinder

    if (r >= radius+cutoff || x[2] <= lo-cutoff || x[2] >= hi+cutoff) return 0;
    if (r < radius && x[2] > lo && x[2] < hi) return 0;

    // z is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    if (r > radius) {
      xp = c1 + del1*radius/r;
      yp = c2 + del2*radius/r;
    } else {
      xp = x[0];
      yp = x[1];
    }
    if (x[2] < lo) zp = lo;
    else if (x[2] > hi) zp = hi;
    else zp = x[2];

    add_contact(0,x,xp,yp,zp);
    if (contact[0].r < cutoff) return 1;
    return 0;
  }
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegCylinder::shape_update()
{
  radius = input->variable->compute_equal(rvar);
  if (radius < 0.0)
    error->one(FLERR,"Variable evaluation in region gave bad value");
  if (axis == 'x') radius *= xscale;
  else radius *= xscale;
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegCylinder::variable_check()
{
  rvar = input->variable->find(rstr);
  if (rvar < 0)
    error->all(FLERR,"Variable name for region cylinder does not exist");
  if (!input->variable->equalstyle(rvar))
    error->all(FLERR,"Variable for region cylinder is invalid style");
}
