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
   Contributing author: Pieter in't Veld (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "region_prism.h"
#include "domain.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

RegPrism::RegPrism(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  options(narg-11,&arg[11]);

  if (strcmp(arg[2],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    xlo = domain->boxxlo;
  } else xlo = xscale*atof(arg[2]);

  if (strcmp(arg[3],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    xhi = domain->boxxhi;
  } else xhi = xscale*atof(arg[3]);

  if (strcmp(arg[4],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    ylo = domain->boxylo;
  } else ylo = yscale*atof(arg[4]);

  if (strcmp(arg[5],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    yhi = domain->boxyhi;
  } else yhi = yscale*atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    zlo = domain->boxzlo;
  } else zlo = zscale*atof(arg[6]);

  if (strcmp(arg[7],"INF") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF when box does not exist");
    zhi = domain->boxzhi;
  } else zhi = zscale*atof(arg[7]);

  yxshift = xscale*atof(arg[8]);
  zxshift = xscale*atof(arg[9]);
  zyshift = yscale*atof(arg[10]);

  // error check
  // prism cannot be 0 thickness in any dim, else inverse blows up

  if (xlo >= xhi || ylo >= yhi || zlo >= zhi)
    error->all("Illegal region prism command");

  // extent of prism
  
  extent_xlo = MIN(xlo,xlo+yxshift);
  extent_xlo = MIN(extent_xlo,extent_xlo+zxshift);
  extent_xhi = MAX(xhi,xhi+yxshift);
  extent_xhi = MAX(extent_xhi,extent_xhi+zxshift);
  extent_ylo = MIN(ylo,ylo+zyshift);
  extent_yhi = MAX(yhi,yhi+zyshift);
  extent_zlo = zlo;
  extent_zhi = zhi;

  // h = transformation matrix from tilt coords (0-1) to box coords (xyz)
  // columns of h are edge vectors of tilted box
  // hinv = transformation matrix from box coords to tilt coords
  // both h and hinv are upper triangular
  //   since 1st edge of prism is along x-axis
  //   and bottom face of prism is in xy plane

  h[0][0] = xhi - xlo;
  h[0][1] = yxshift;
  h[1][1] = yhi - ylo;
  h[0][2] = zxshift;
  h[1][2] = zyshift;
  h[2][2] = zhi - zlo;

  hinv[0][0] = 1.0/h[0][0];
  hinv[0][1] = -h[0][1] / (h[0][0]*h[1][1]);
  hinv[1][1] = 1.0/h[1][1];
  hinv[0][2] = (h[0][1]*h[1][2] - h[0][2]*h[1][1]) / (h[0][0]*h[1][1]*h[2][2]);
  hinv[1][2] = -h[1][2] / (h[1][1]*h[2][2]);
  hinv[2][2] = 1.0/h[2][2];
}

/* ----------------------------------------------------------------------
   check xyz against prism
   abc = Hinv * (xyz - xyzlo)
   abc = tilt coords (0-1)
   Hinv = transformation matrix from box coords to tilt coords
   xyz = box coords
   xyzlo = lower-left corner of prism
------------------------------------------------------------------------- */

int RegPrism::match(double x, double y, double z)
{
  double a = hinv[0][0]*(x-xlo) + hinv[0][1]*(y-ylo) + hinv[0][2]*(z-zlo);
  double b = hinv[1][1]*(y-ylo) + hinv[1][2]*(z-zlo);
  double c = hinv[2][2]*(z-zlo);

  int inside;
  if (a >= 0.0 && a <= 1.0 && b >= 0.0 && b <= 1.0 && c >= 0.0 && c <= 1.0)
    inside = 1;
  else inside = 0;

  return !(inside ^ interior);         // 1 if same, 0 if different
}
