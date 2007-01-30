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
#include "region_sphere.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegSphere::RegSphere(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  options(narg-6,&arg[6]);

  xc = xscale*atof(arg[2]);
  yc = yscale*atof(arg[3]);
  zc = zscale*atof(arg[4]);
  radius = xscale*atof(arg[5]);

  // error check

  if (radius < 0.0) error->all("Illegal region sphere command");

  // extent of sphere

  extent_xlo = xc - radius;
  extent_xhi = xc + radius;
  extent_ylo = yc - radius;
  extent_yhi = yc + radius;
  extent_zlo = zc - radius;
  extent_zhi = zc + radius;
}

/* ---------------------------------------------------------------------- */

int RegSphere::match(double x, double y, double z)
{
  double delx = x - xc;
  double dely = y - yc;
  double delz = z - zc;
  double dist = sqrt(delx*delx + dely*dely + delz*delz);

  int inside;
  if (dist <= radius) inside = 1;
  else inside = 0;

  return !(inside ^ interior);         // 1 if same, 0 if different
}
