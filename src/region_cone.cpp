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
  radiuslo = xscale*atof(arg[5]);
  radiushi = xscale*atof(arg[6]);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
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
      error->all("Cannot use region INF or EDGE when box does not exist");
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

  if (radiuslo < 0.0) error->all("Illegal radius in region cone command");
  if (radiushi < 0.0) error->all("Illegal radius in region cone command");
  if (hi == lo) error->all("Illegal cone length in region cone command");

  // extent of cone

  double largestradius;
  if (radiuslo > radiushi) largestradius = radiuslo;
  else largestradius = radiushi;

  if (axis == 'x') {
    extent_xlo = lo;
    extent_xhi = hi;
    extent_ylo = c1 - largestradius;
    extent_yhi = c1 + largestradius;
    extent_zlo = c2 - largestradius;
    extent_zhi = c2 + largestradius;
  }
  if (axis == 'y') {
    extent_xlo = c1 - largestradius;
    extent_xhi = c1 + largestradius;
    extent_ylo = lo;
    extent_yhi = hi;
    extent_zlo = c2 - largestradius;
    extent_zhi = c2 + largestradius;
  }
  if (axis == 'z') {
    extent_xlo = c1 - largestradius;
    extent_xhi = c1 + largestradius;
    extent_ylo = c2 - largestradius;
    extent_yhi = c2 + largestradius;
    extent_zlo = lo;
    extent_zhi = hi;
  }
}

/* ---------------------------------------------------------------------- */

int RegCone::match(double x, double y, double z)
{
  double del1,del2,dist;
  int inside;
  double currentradius;

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

  return !(inside ^ interior);         // 1 if same, 0 if different
}
