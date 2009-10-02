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
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegCylinder::RegCylinder(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"x") && strcmp(arg[2],"y") && strcmp(arg[2],"z")) 
    error->all("Illegal region cylinder command");
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
  radius = xscale*atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
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
      error->all("Cannot use region INF or EDGE when box does not exist");
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

  if (radius < 0.0) error->all("Illegal region cylinder command");

  // extent of cylinder

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
}

/* ---------------------------------------------------------------------- */

int RegCylinder::match(double x, double y, double z)
{
  double del1,del2,dist;
  int inside;

  if (axis == 'x') {
    del1 = y - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && x >= lo && x <= hi) inside = 1;
    else inside = 0;
  }
  if (axis == 'y') {
    del1 = x - c1;
    del2 = z - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && y >= lo && y <= hi) inside = 1;
    else inside = 0;
  }
  if (axis == 'z') {
    del1 = x - c1;
    del2 = y - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && z >= lo && z <= hi) inside = 1;
    else inside = 0;
  }

  return !(inside ^ interior);         // 1 if same, 0 if different
}
