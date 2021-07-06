// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "region_block.h"

#include "domain.h"
#include "error.h"
#include "math_extra.h"

#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegBlock::RegBlock(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else if (domain->triclinic == 0) xlo = domain->boxlo[0];
    else xlo = domain->boxlo_bound[0];
  } else xlo = xscale*utils::numeric(FLERR,arg[2],false,lmp);

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else if (domain->triclinic == 0) xhi = domain->boxhi[0];
    else xhi = domain->boxhi_bound[0];
  } else xhi = xscale*utils::numeric(FLERR,arg[3],false,lmp);

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else if (domain->triclinic == 0) ylo = domain->boxlo[1];
    else ylo = domain->boxlo_bound[1];
  } else ylo = yscale*utils::numeric(FLERR,arg[4],false,lmp);

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else if (domain->triclinic == 0) yhi = domain->boxhi[1];
    else yhi = domain->boxhi_bound[1];
  } else yhi = yscale*utils::numeric(FLERR,arg[5],false,lmp);

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else if (domain->triclinic == 0) zlo = domain->boxlo[2];
    else zlo = domain->boxlo_bound[2];
  } else zlo = zscale*utils::numeric(FLERR,arg[6],false,lmp);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else if (domain->triclinic == 0) zhi = domain->boxhi[2];
    else zhi = domain->boxhi_bound[2];
  } else zhi = zscale*utils::numeric(FLERR,arg[7],false,lmp);

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"Illegal region block command");

  // extent of block

  if (interior) {
    bboxflag = 1;
    extent_xlo = xlo;
    extent_xhi = xhi;
    extent_ylo = ylo;
    extent_yhi = yhi;
    extent_zlo = zlo;
    extent_zhi = zhi;
  } else bboxflag = 0;

  // particle could be close to all 6 planes
  // particle can only touch 3 planes

  cmax = 6;
  contact = new Contact[cmax];
  if (interior) tmax = 3;
  else tmax = 1;

  // open face data structs

  face[0][0] = -1.0;
  face[0][1] = 0.0;
  face[0][2] = 0.0;
  face[1][0] = 1.0;
  face[1][1] = 0.0;
  face[1][2] = 0.0;
  face[2][0] = 0.0;
  face[2][1] = -1.0;
  face[2][2] = 0.0;
  face[3][0] = 0.0;
  face[3][1] = 1.0;
  face[3][2] = 0.0;
  face[4][0] = 0.0;
  face[4][1] = 0.0;
  face[4][2] = -1.0;
  face[5][0] = 0.0;
  face[5][1] = 0.0;
  face[5][2] = 1.0;

  // face[0]

  corners[0][0][0] = xlo;
  corners[0][0][1] = ylo;
  corners[0][0][2] = zlo;
  corners[0][1][0] = xlo;
  corners[0][1][1] = ylo;
  corners[0][1][2] = zhi;
  corners[0][2][0] = xlo;
  corners[0][2][1] = yhi;
  corners[0][2][2] = zhi;
  corners[0][3][0] = xlo;
  corners[0][3][1] = yhi;
  corners[0][3][2] = zlo;

  // face[1]

  corners[1][0][0] = xhi;
  corners[1][0][1] = ylo;
  corners[1][0][2] = zlo;
  corners[1][1][0] = xhi;
  corners[1][1][1] = ylo;
  corners[1][1][2] = zhi;
  corners[1][2][0] = xhi;
  corners[1][2][1] = yhi;
  corners[1][2][2] = zhi;
  corners[1][3][0] = xhi;
  corners[1][3][1] = yhi;
  corners[1][3][2] = zlo;

  // face[2]

  MathExtra::copy3(corners[0][0], corners[2][0]);
  MathExtra::copy3(corners[1][0], corners[2][1]);
  MathExtra::copy3(corners[1][1], corners[2][2]);
  MathExtra::copy3(corners[0][1], corners[2][3]);

  // face[3]

  MathExtra::copy3(corners[0][3], corners[3][0]);
  MathExtra::copy3(corners[0][2], corners[3][1]);
  MathExtra::copy3(corners[1][2], corners[3][2]);
  MathExtra::copy3(corners[1][3], corners[3][3]);

  // face[4]

  MathExtra::copy3(corners[0][0], corners[4][0]);
  MathExtra::copy3(corners[0][3], corners[4][1]);
  MathExtra::copy3(corners[1][3], corners[4][2]);
  MathExtra::copy3(corners[1][0], corners[4][3]);

  // face[5]

  MathExtra::copy3(corners[0][1], corners[5][0]);
  MathExtra::copy3(corners[1][1], corners[5][1]);
  MathExtra::copy3(corners[1][2], corners[5][2]);
  MathExtra::copy3(corners[0][2], corners[5][3]);
}

/* ---------------------------------------------------------------------- */

RegBlock::~RegBlock()
{
  if (copymode) return;

  delete [] contact;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegBlock::inside(double x, double y, double z)
{
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of block
   can be one contact for each of 6 faces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

int RegBlock::surface_interior(double *x, double cutoff)
{
  double delta;

  // x is exterior to block

  if (x[0] < xlo || x[0] > xhi || x[1] < ylo || x[1] > yhi ||
      x[2] < zlo || x[2] > zhi) return 0;

  // x is interior to block or on its surface

  int n = 0;

  delta = x[0] - xlo;
  if (delta < cutoff && !open_faces[0]) {
    contact[n].r = delta;
    contact[n].delx = delta;
    contact[n].dely = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 0;
    n++;
  }
  delta = xhi - x[0];
  if (delta < cutoff && !open_faces[1]) {
    contact[n].r = delta;
    contact[n].delx = -delta;
    contact[n].dely = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 1;
    n++;
  }

  delta = x[1] - ylo;
  if (delta < cutoff && !open_faces[2]) {
    contact[n].r = delta;
    contact[n].dely = delta;
    contact[n].delx = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 2;
    n++;
  }
  delta = yhi - x[1];
  if (delta < cutoff && !open_faces[3]) {
    contact[n].r = delta;
    contact[n].dely = -delta;
    contact[n].delx = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 3;
    n++;
  }

  delta = x[2] - zlo;
  if (delta < cutoff && !open_faces[4]) {
    contact[n].r = delta;
    contact[n].delz = delta;
    contact[n].delx = contact[n].dely = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 4;
    n++;
  }
  delta = zhi - x[2];
  if (delta < cutoff && !open_faces[5]) {
    contact[n].r = delta;
    contact[n].delz = -delta;
    contact[n].delx = contact[n].dely = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 5;
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of block
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

int RegBlock::surface_exterior(double *x, double cutoff)
{
  double xp,yp,zp;
  double xc,yc,zc,dist,mindist;

  // x is far enough from block that there is no contact
  // x is interior to block

  if (x[0] <= xlo-cutoff || x[0] >= xhi+cutoff ||
      x[1] <= ylo-cutoff || x[1] >= yhi+cutoff ||
      x[2] <= zlo-cutoff || x[2] >= zhi+cutoff) return 0;
  if (x[0] > xlo && x[0] < xhi && x[1] > ylo && x[1] < yhi &&
      x[2] > zlo && x[2] < zhi) return 0;

  // x is exterior to block or on its surface
  // xp,yp,zp = point on surface of block that x is closest to
  //            could be edge or corner pt of block
  // do not add contact point if r >= cutoff

  if (!openflag) {
    if (x[0] < xlo) xp = xlo;
    else if (x[0] > xhi) xp = xhi;
    else xp = x[0];
    if (x[1] < ylo) yp = ylo;
    else if (x[1] > yhi) yp = yhi;
    else yp = x[1];
    if (x[2] < zlo) zp = zlo;
    else if (x[2] > zhi) zp = zhi;
    else zp = x[2];
  } else {
    mindist = BIG;
    for (int i = 0; i < 6; i++) {
      if (open_faces[i]) continue;
      dist = find_closest_point(i,x,xc,yc,zc);
      if (dist < mindist) {
        xp = xc;
        yp = yc;
        zp = zc;
        mindist = dist;
      }
    }
  }

  add_contact(0,x,xp,yp,zp);
  contact[0].iwall = 0;
  if (contact[0].r < cutoff) return 1;
  return 0;
}

/*------------------------------------------------------------------------
  return distance to closest point on surface I of block region
  store closest point in xc,yc,zc
--------------------------------------------------------------------------*/

double RegBlock::find_closest_point(int i, double *x,
                                    double &xc, double &yc, double &zc)
{
  double dot,d2,d2min;
  double xr[3],xproj[3],p[3];

  xr[0] = x[0] - corners[i][0][0];
  xr[1] = x[1] - corners[i][0][1];
  xr[2] = x[2] - corners[i][0][2];
  dot = face[i][0]*xr[0] + face[i][1]*xr[1] + face[i][2]*xr[2];
  xproj[0] = xr[0] - dot*face[i][0];
  xproj[1] = xr[1] - dot*face[i][1];
  xproj[2] = xr[2] - dot*face[i][2];

  d2min = BIG;

  // check if point projects inside of face

  if (inside_face(xproj, i)) {
    d2 = d2min = dot*dot;
    xc = xproj[0] + corners[i][0][0];
    yc = xproj[1] + corners[i][0][1];
    zc = xproj[2] + corners[i][0][2];

 // check each edge

  } else {
    point_on_line_segment(corners[i][0],corners[i][1],x,p);
    d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
      (p[2]-x[2])*(p[2]-x[2]);
    if (d2 < d2min) {
      d2min = d2;
      xc = p[0];
      yc = p[1];
      zc = p[2];
    }

    point_on_line_segment(corners[i][1],corners[i][2],x,p);
    d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
      (p[2]-x[2])*(p[2]-x[2]);
    if (d2 < d2min) {
      d2min = d2;
      xc = p[0];
      yc = p[1];
      zc = p[2];
    }

    point_on_line_segment(corners[i][2],corners[i][3],x,p);
    d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
      (p[2]-x[2])*(p[2]-x[2]);
    if (d2 < d2min) {
      d2min = d2;
      xc = p[0];
      yc = p[1];
      zc = p[2];
    }

    point_on_line_segment(corners[i][3],corners[i][0],x,p);
    d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
      (p[2]-x[2])*(p[2]-x[2]);
    if (d2 < d2min) {
      d2min = d2;
      xc = p[0];
      yc = p[1];
      zc = p[2];
    }
  }

  return d2min;
}

/*------------------------------------------------------------------------
  determine if projected point is inside given face of the block
--------------------------------------------------------------------------*/

int RegBlock::inside_face(double *xproj, int iface)
{
  if (iface < 2) {
    if (xproj[1] > 0 && (xproj[1] < yhi-ylo) &&
        xproj[2] > 0 && (xproj[2] < zhi-zlo)) return 1;
  } else if (iface < 4) {
    if (xproj[0] > 0 && (xproj[0] < (xhi-xlo)) &&
        xproj[2] > 0 && (xproj[2] < (zhi-zlo))) return 1;
  } else {
    if (xproj[0] > 0 && xproj[0] < (xhi-xlo) &&
        xproj[1] > 0 && xproj[1] < (yhi-ylo)) return 1;
  }

  return 0;
}
