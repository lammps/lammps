/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "region.h"
#include "domain.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Region::Region(int narg, char **arg)
{
  // region ID

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  // set option defaults

  interior = 1;
  scaleflag = 1;

  // read style and options from end of input line

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  if (strcmp(style,"block") == 0) options(narg-8,&arg[8]);
  else if (strcmp(style,"sphere") == 0) options(narg-6,&arg[6]);
  else if (strcmp(arg[1],"cylinder") == 0) options(narg-8,&arg[8]);
  else if (strcmp(style,"prism") == 0) options(narg-11,&arg[11]);
  else if (strcmp(arg[1],"union") == 0) {
    if (narg < 5) error->all("Illegal region command");
    n = atoi(arg[2]);
    options(narg-(n+3),&arg[n+3]);
  } else if (strcmp(arg[1],"intersect") == 0) {
    if (narg < 5) error->all("Illegal region command");
    n = atoi(arg[2]);
    options(narg-(n+3),&arg[n+3]);
  } else error->all("Illegal region command");

  // setup scaling

  if (scaleflag && strcmp(domain->lattice_style,"none") == 0)
    error->all("Use of region with undefined lattice");

  if (scaleflag) {
    xscale = domain->xlattice;
    yscale = domain->ylattice;
    zscale = domain->zlattice;
  }
  else xscale = yscale = zscale = 1.0;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  delete [] id;
  delete [] style;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal region command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal region command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all("Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all("Illegal region command");
      iarg += 2;
    } else error->all("Illegal region command");
  }
}
