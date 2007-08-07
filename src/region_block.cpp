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

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegBlock::RegBlock(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else if (domain->triclinic == 0) xlo = domain->boxlo[0];
    else xlo = domain->boxlo_bound[0];
  } else xlo = xscale*atof(arg[2]);

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else if (domain->triclinic == 0) xhi = domain->boxhi[0];
    else xhi = domain->boxhi_bound[0];
  } else xhi = xscale*atof(arg[3]);

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else if (domain->triclinic == 0) ylo = domain->boxlo[1];
    else ylo = domain->boxlo_bound[1];
  } else ylo = yscale*atof(arg[4]);

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else if (domain->triclinic == 0) yhi = domain->boxhi[1];
    else yhi = domain->boxhi_bound[1];
  } else yhi = yscale*atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else if (domain->triclinic == 0) zlo = domain->boxlo[2];
    else zlo = domain->boxlo_bound[2];
  } else zlo = zscale*atof(arg[6]);

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0) 
      error->all("Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else if (domain->triclinic == 0) zhi = domain->boxhi[2];
    else zhi = domain->boxhi_bound[2];
  } else zhi = zscale*atof(arg[7]);

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all("Illegal region block command");

  // extent of block
  
  extent_xlo = xlo;
  extent_xhi = xhi;
  extent_ylo = ylo;
  extent_yhi = yhi;
  extent_zlo = zlo;
  extent_zhi = zhi;
}

/* ---------------------------------------------------------------------- */

int RegBlock::match(double x, double y, double z)
{
  int inside;
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    inside = 1;
  else inside = 0;

  return !(inside ^ interior);         // 1 if same, 0 if different
}
