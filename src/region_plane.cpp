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
#include "region_plane.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegPlane::RegPlane(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);

  xp = xscale*force->numeric(FLERR,arg[2]);
  yp = yscale*force->numeric(FLERR,arg[3]);
  zp = zscale*force->numeric(FLERR,arg[4]);
  normal[0] = xscale*force->numeric(FLERR,arg[5]);
  normal[1] = yscale*force->numeric(FLERR,arg[6]);
  normal[2] = zscale*force->numeric(FLERR,arg[7]);

  // enforce unit normal

  double rsq = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
  if (rsq == 0.0) error->all(FLERR,"Illegal region plane command");
  normal[0] /= sqrt(rsq);
  normal[1] /= sqrt(rsq);
  normal[2] /= sqrt(rsq);

  // plane has no bounding box

  bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegPlane::~RegPlane()
{
  delete [] contact;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is on normal side of plane or on plane
   inside = 0 if x,y,z is on non-normal side of plane and not on plane
   x,y,z is inside if (x-xp) dot normal >= 0
------------------------------------------------------------------------- */

int RegPlane::inside(double x, double y, double z)
{
  double dot = (x-xp)*normal[0] + (y-yp)*normal[1] + (z-zp)*normal[2];

  if (dot >= 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from normal side of plane
   no contact if on other side (possible if called from union/intersect)
   delxyz = vector from nearest projected point on plane to x
------------------------------------------------------------------------- */

int RegPlane::surface_interior(double *x, double cutoff)
{
  double dot = (x[0]-xp)*normal[0] + (x[1]-yp)*normal[1] + (x[2]-zp)*normal[2];
  if (dot < cutoff && dot >= 0.0) {
    contact[0].r = dot;
    contact[0].delx = dot*normal[0];
    contact[0].dely = dot*normal[1];
    contact[0].delz = dot*normal[2];
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from non-normal side of plane
   no contact if on other side (possible if called from union/intersect)
   delxyz = vector from nearest projected point on plane to x
------------------------------------------------------------------------- */

int RegPlane::surface_exterior(double *x, double cutoff)
{
  double dot = (x[0]-xp)*normal[0] + (x[1]-yp)*normal[1] + (x[2]-zp)*normal[2];
  dot = -dot;
  if (dot < cutoff && dot >= 0.0) {
    contact[0].r = dot;
    contact[0].delx = -dot*normal[0];
    contact[0].dely = -dot*normal[1];
    contact[0].delz = -dot*normal[2];
    return 1;
  }
  return 0;
}
