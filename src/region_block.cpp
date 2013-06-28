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

#include "stdlib.h"
#include "string.h"
#include "region_block.h"
#include "domain.h"
#include "error.h"
#include "force.h"

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
  } else xlo = xscale*force->numeric(FLERR,arg[2]);

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else if (domain->triclinic == 0) xhi = domain->boxhi[0];
    else xhi = domain->boxhi_bound[0];
  } else xhi = xscale*force->numeric(FLERR,arg[3]);

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else if (domain->triclinic == 0) ylo = domain->boxlo[1];
    else ylo = domain->boxlo_bound[1];
  } else ylo = yscale*force->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else if (domain->triclinic == 0) yhi = domain->boxhi[1];
    else yhi = domain->boxhi_bound[1];
  } else yhi = yscale*force->numeric(FLERR,arg[5]);

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else if (domain->triclinic == 0) zlo = domain->boxlo[2];
    else zlo = domain->boxlo_bound[2];
  } else zlo = zscale*force->numeric(FLERR,arg[6]);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else if (domain->triclinic == 0) zhi = domain->boxhi[2];
    else zhi = domain->boxhi_bound[2];
  } else zhi = zscale*force->numeric(FLERR,arg[7]);

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

  cmax = 6;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegBlock::~RegBlock()
{
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
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delx = delta;
    contact[n].dely = contact[n].delz = 0.0;
    n++;
  }
  delta = xhi - x[0];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delx = -delta;
    contact[n].dely = contact[n].delz = 0.0;
    n++;
  }

  delta = x[1] - ylo;
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].dely = delta;
    contact[n].delx = contact[n].delz = 0.0;
    n++;
  }
  delta = yhi - x[1];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].dely = -delta;
    contact[n].delx = contact[n].delz = 0.0;
    n++;
  }

  delta = x[2] - zlo;
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delz = delta;
    contact[n].delx = contact[n].dely = 0.0;
    n++;
  }
  delta = zhi - x[2];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delz = -delta;
    contact[n].delx = contact[n].dely = 0.0;
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

  if (x[0] < xlo) xp = xlo;
  else if (x[0] > xhi) xp = xhi;
  else xp = x[0];
  if (x[1] < ylo) yp = ylo;
  else if (x[1] > yhi) yp = yhi;
  else yp = x[1];
  if (x[2] < zlo) zp = zlo;
  else if (x[2] > zhi) zp = zhi;
  else zp = x[2];

  add_contact(0,x,xp,yp,zp);
  if (contact[0].r < cutoff) return 1;
  return 0;
}
