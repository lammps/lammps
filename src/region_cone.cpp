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
   Contributing author: Pim Schravendijk
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "region_cone.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegCone::RegCone(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-9,&arg[9]);

  if (strcmp(arg[2],"x") && strcmp(arg[2],"y") && strcmp(arg[2],"z"))
    error->all(FLERR,"Illegal region cylinder command");
  axis = arg[2][0];

  if (axis == 'x') {
    c1 = yscale*atof(arg[3]);
    c2 = zscale*atof(arg[4]);
    radiuslo = yscale*atof(arg[5]);
    radiushi = yscale*atof(arg[6]);
  } else if (axis == 'y') {
    c1 = xscale*atof(arg[3]);
    c2 = zscale*atof(arg[4]);
    radiuslo = xscale*atof(arg[5]);
    radiushi = xscale*atof(arg[6]);
  } else if (axis == 'z') {
    c1 = xscale*atof(arg[3]);
    c2 = yscale*atof(arg[4]);
    radiuslo = xscale*atof(arg[5]);
    radiushi = xscale*atof(arg[6]);
  }

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[7],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[0];
      else lo = domain->boxlo_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[7],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[1];
      else lo = domain->boxlo_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[7],"INF") == 0) lo = -BIG;
      else if (domain->triclinic == 0) lo = domain->boxlo[2];
      else lo = domain->boxlo_bound[2];
    }
  } else {
    if (axis == 'x') lo = xscale*atof(arg[7]);
    if (axis == 'y') lo = yscale*atof(arg[7]);
    if (axis == 'z') lo = zscale*atof(arg[7]);
  }

  if (strcmp(arg[8],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[8],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[0];
      else hi = domain->boxhi_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[8],"INF") == 0) hi = BIG;
      if (domain->triclinic == 0) hi = domain->boxhi[1];
      else hi = domain->boxhi_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[8],"INF") == 0) hi = BIG;
      else if (domain->triclinic == 0) hi = domain->boxhi[2];
      else hi = domain->boxhi_bound[2];
    }
  } else {
    if (axis == 'x') hi = xscale*atof(arg[8]);
    if (axis == 'y') hi = yscale*atof(arg[8]);
    if (axis == 'z') hi = zscale*atof(arg[8]);
  }

  // error check

  if (radiuslo < 0.0) error->all(FLERR,"Illegal radius in region cone command");
  if (radiushi < 0.0) error->all(FLERR,"Illegal radius in region cone command");
  if (radiuslo == 0.0 && radiushi == 0.0)
    error->all(FLERR,"Illegal radius in region cone command");
  if (hi == lo) error->all(FLERR,"Illegal cone length in region cone command");

  // extent of cone

  if (radiuslo > radiushi) maxradius = radiuslo;
  else maxradius = radiushi;

  if (interior) {
    bboxflag = 1;

    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 - maxradius;
      extent_yhi = c1 + maxradius;
      extent_zlo = c2 - maxradius;
      extent_zhi = c2 + maxradius;
    }
    if (axis == 'y') {
      extent_xlo = c1 - maxradius;
      extent_xhi = c1 + maxradius;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 - maxradius;
      extent_zhi = c2 + maxradius;
    }
    if (axis == 'z') {
      extent_xlo = c1 - maxradius;
      extent_xhi = c1 + maxradius;
      extent_ylo = c2 - maxradius;
      extent_yhi = c2 + maxradius;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else bboxflag = 0;

  // particle could be contact cone surface and 2 ends

  cmax = 3;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegCone::~RegCone()
{
  delete [] contact;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCone::inside(double x, double y, double z)
{
  double del1,del2,dist;
  double currentradius;
  int inside;

  if (axis == 'x') {
    del1 = y - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x-lo)*(radiushi-radiuslo)/(hi-lo);
    if (dist <= currentradius && x >= lo && x <= hi) inside = 1;
    else inside = 0;
  }
  if (axis == 'y') {
    del1 = x - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (y-lo)*(radiushi-radiuslo)/(hi-lo);
    if (dist <= currentradius && y >= lo && y <= hi) inside = 1;
    else inside = 0;
  }
  if (axis == 'z') {
    del1 = x - c1;
    del2 = y - c2;
    dist = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (z-lo)*(radiushi-radiuslo)/(hi-lo);
    if (dist <= currentradius && z >= lo && z <= hi) inside = 1;
    else inside = 0;
  }

  return inside;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of cone
   can be one contact for each of 3 cone surfaces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on cone to x
   special case: no contact with curved surf if x is on center axis
------------------------------------------------------------------------- */

int RegCone::surface_interior(double *x, double cutoff)
{
  double del1,del2,r,currentradius,delx,dely,delz,dist,delta;
  double surflo[3],surfhi[3],xs[3];

  int n = 0;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x[0]-lo)*(radiushi-radiuslo)/(hi-lo);

    // x is exterior to cone

    if (r > currentradius || x[0] < lo || x[0] > hi) return 0;

    // x is interior to cone or on its surface
    // surflo = pt on outer circle of bottom end plane, same dir as x vs axis
    // surfhi = pt on outer circle of top end plane, same dir as x vs axis

    if (r > 0.0) {
      surflo[0] = lo;
      surflo[1] = c1 + del1*radiuslo/r;
      surflo[2] = c2 + del2*radiuslo/r;
      surfhi[0] = hi;
      surfhi[1] = c1 + del1*radiushi/r;
      surfhi[2] = c2 + del2*radiushi/r;
      point_on_line_segment(surflo,surfhi,x,xs);
      delx = x[0] - xs[0];
      dely = x[1] - xs[1];
      delz = x[2] - xs[2];
      dist = sqrt(delx*delx + dely*dely + delz*delz);
      if (dist < cutoff) {
        contact[n].r = dist;
        contact[n].delx = delx;
        contact[n].dely = dely;
        contact[n].delz = delz;
        n++;
      }
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
    currentradius = radiuslo + (x[1]-lo)*(radiushi-radiuslo)/(hi-lo);

    // y is exterior to cone

    if (r > currentradius || x[1] < lo || x[1] > hi) return 0;

    // y is interior to cone or on its surface
    // surflo = pt on outer circle of bottom end plane, same dir as y vs axis
    // surfhi = pt on outer circle of top end plane, same dir as y vs axis

    if (r > 0.0) {
      surflo[0] = c1 + del1*radiuslo/r;
      surflo[1] = lo;
      surflo[2] = c2 + del2*radiuslo/r;
      surfhi[0] = c1 + del1*radiushi/r;
      surfhi[1] = hi;
      surfhi[2] = c2 + del2*radiushi/r;
      point_on_line_segment(surflo,surfhi,x,xs);
      delx = x[0] - xs[0];
      dely = x[1] - xs[1];
      delz = x[2] - xs[2];
      dist = sqrt(delx*delx + dely*dely + delz*delz);
      if (dist < cutoff) {
        contact[n].r = dist;
        contact[n].delx = delx;
        contact[n].dely = dely;
        contact[n].delz = delz;
        n++;
      }
    }

    delta = x[1] - lo;
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delz = delta;
      contact[n].delx = contact[n].dely = 0.0;
      n++;
    }
    delta = hi - x[1];
    if (delta < cutoff) {
      contact[n].r = delta;
      contact[n].delz = -delta;
      contact[n].delx = contact[n].dely = 0.0;
      n++;
    }

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x[2]-lo)*(radiushi-radiuslo)/(hi-lo);

    // z is exterior to cone

    if (r > currentradius || x[2] < lo || x[2] > hi) return 0;

    // z is interior to cone or on its surface
    // surflo = pt on outer circle of bottom end plane, same dir as z vs axis
    // surfhi = pt on outer circle of top end plane, same dir as z vs axis

    if (r > 0.0) {
      surflo[0] = c1 + del1*radiuslo/r;
      surflo[1] = c2 + del2*radiuslo/r;
      surflo[2] = lo;
      surfhi[0] = c1 + del1*radiushi/r;
      surfhi[1] = c2 + del2*radiushi/r;
      surfhi[2] = hi;
      point_on_line_segment(surflo,surfhi,x,xs);
      delx = x[0] - xs[0];
      dely = x[1] - xs[1];
      delz = x[2] - xs[2];
      dist = sqrt(delx*delx + dely*dely + delz*delz);
      if (dist < cutoff) {
        contact[n].r = dist;
        contact[n].delx = delx;
        contact[n].dely = dely;
        contact[n].delz = delz;
        n++;
      }
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
   one contact if 0 <= x < cutoff from outer surface of cone
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on cone to x
------------------------------------------------------------------------- */

int RegCone::surface_exterior(double *x, double cutoff)
{
  double del1,del2,r,currentradius,distsq;
  double corner1[3],corner2[3],corner3[3],corner4[3],xp[3],nearest[3];

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x[0]-lo)*(radiushi-radiuslo)/(hi-lo);

    // x is far enough from cone that there is no contact
    // x is interior to cone

    if (r >= maxradius+cutoff ||
        x[0] <= lo-cutoff || x[0] >= hi+cutoff) return 0;
    if (r < currentradius && x[0] > lo && x[0] < hi) return 0;

    // x is exterior to cone or on its surface
    // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of x
    // project x to 3 line segments in half trapezoid (4th is axis of cone)
    // nearest = point on surface of cone that x is closest to
    //           could be edge of cone
    // do not add contact point if r >= cutoff

    corner1[0] = lo;
    corner1[1] = c1 + del1*radiuslo/r;
    corner1[2] = c2 + del2*radiuslo/r;
    corner2[0] = hi;
    corner2[1] = c1 + del1*radiushi/r;
    corner2[2] = c2 + del2*radiushi/r;
    corner3[0] = lo;
    corner3[1] = c1;
    corner3[2] = c2;
    corner4[0] = hi;
    corner4[1] = c1;
    corner4[2] = c2;

    distsq = BIG;
    point_on_line_segment(corner1,corner2,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner1,corner3,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner2,corner4,x,xp);
    distsq = closest(x,xp,nearest,distsq);

    add_contact(0,x,nearest[0],nearest[1],nearest[2]);
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x[1]-lo)*(radiushi-radiuslo)/(hi-lo);

    // y is far enough from cone that there is no contact
    // y is interior to cone

    if (r >= maxradius+cutoff ||
        x[1] <= lo-cutoff || x[1] >= hi+cutoff) return 0;
    if (r < currentradius && x[1] > lo && x[1] < hi) return 0;

    // y is exterior to cone or on its surface
    // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of y
    // project x to 3 line segments in half trapezoid (4th is axis of cone)
    // nearest = point on surface of cone that y is closest to
    //           could be edge of cone
    // do not add contact point if r >= cutoff

    corner1[0] = c1 + del1*radiuslo/r;
    corner1[1] = lo;
    corner1[2] = c2 + del2*radiuslo/r;
    corner2[0] = c1 + del1*radiushi/r;
    corner2[1] = hi;
    corner2[2] = c2 + del2*radiushi/r;
    corner3[0] = c1;
    corner3[1] = lo;
    corner3[2] = c2;
    corner4[0] = c1;
    corner4[1] = hi;
    corner4[2] = c2;

    distsq = BIG;
    point_on_line_segment(corner1,corner2,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner1,corner3,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner2,corner4,x,xp);
    distsq = closest(x,xp,nearest,distsq);

    add_contact(0,x,nearest[0],nearest[1],nearest[2]);
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1*del1 + del2*del2);
    currentradius = radiuslo + (x[2]-lo)*(radiushi-radiuslo)/(hi-lo);

    // z is far enough from cone that there is no contact
    // z is interior to cone

    if (r >= maxradius+cutoff ||
        x[2] <= lo-cutoff || x[2] >= hi+cutoff) return 0;
    if (r < currentradius && x[2] > lo && x[2] < hi) return 0;

    // z is exterior to cone or on its surface
    // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of z
    // project x to 3 line segments in half trapezoid (4th is axis of cone)
    // nearest = point on surface of cone that z is closest to
    //           could be edge of cone
    // do not add contact point if r >= cutoff

    corner1[0] = c1 + del1*radiuslo/r;
    corner1[1] = c2 + del2*radiuslo/r;
    corner1[2] = lo;
    corner2[0] = c1 + del1*radiushi/r;
    corner2[1] = c2 + del2*radiushi/r;
    corner2[2] = hi;
    corner3[0] = c1;
    corner3[1] = c2;
    corner3[2] = lo;
    corner4[0] = c1;
    corner4[1] = c2;
    corner4[2] = hi;

    distsq = BIG;
    point_on_line_segment(corner1,corner2,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner1,corner3,x,xp);
    distsq = closest(x,xp,nearest,distsq);
    point_on_line_segment(corner2,corner4,x,xp);
    distsq = closest(x,xp,nearest,distsq);

    add_contact(0,x,nearest[0],nearest[1],nearest[2]);
    if (contact[0].r < cutoff) return 1;
    return 0;
  }
}

/* ----------------------------------------------------------------------
   find nearest point to C on line segment A,B and return it as D
   project (C-A) onto (B-A)
   t = length of that projection, normalized by length of (B-A)
   t <= 0, C is closest to A
   t >= 1, C is closest to B
   else closest point is between A and B
------------------------------------------------------------------------- */

void RegCone::point_on_line_segment(double *a, double *b,
                                    double *c, double *d)
{
  double ba[3],ca[3];

  subtract(a,b,ba);
  subtract(a,c,ca);
  double t = dotproduct(ca,ba) / dotproduct(ba,ba);

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

/* ---------------------------------------------------------------------- */

double RegCone::closest(double *x, double *near, double *nearest, double dsq)
{
  double delx = x[0] - near[0];
  double dely = x[1] - near[1];
  double delz = x[2] - near[2];
  double rsq = delx*delx + dely*dely + delz*delz;
  if (rsq >= dsq) return dsq;

  nearest[0] = near[0];
  nearest[1] = near[1];
  nearest[2] = near[2];
  return rsq;
}

/* ----------------------------------------------------------------------
   v3 = v2 - v1
------------------------------------------------------------------------- */

void RegCone::subtract(double *v1, double *v2, double *v3)
{
  v3[0] = v2[0] - v1[0];
  v3[1] = v2[1] - v1[1];
  v3[2] = v2[2] - v1[2];
}

/* ----------------------------------------------------------------------
   return dotproduct = v1 dot v2
------------------------------------------------------------------------- */

double RegCone::dotproduct(double *v1, double *v2)
{
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}
