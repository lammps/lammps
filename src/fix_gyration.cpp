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
#include "fix_gyration.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixGyration::FixGyration(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal fix gyration command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix gyration command");
  first = 1;

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix gyration file %s",arg[4]);
      error->one(str);
    }
  }

  if (me == 0) {
    fprintf(fp,"# Radius-of-gyration for group %s\n",group->names[igroup]);
    fprintf(fp,"# TimeStep Rg\n");
  }
}

/* ---------------------------------------------------------------------- */

FixGyration::~FixGyration()
{
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixGyration::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGyration::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixGyration::setup(int vflag)
{
  if (first) end_of_step();
  first = 0;
}

/* ---------------------------------------------------------------------- */

void FixGyration::end_of_step()
{
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);
  double rg = group->gyration(igroup,masstotal,xcm);

  if (me == 0) {
    fprintf(fp,"%d %g\n",update->ntimestep,rg);
    fflush(fp);
  }
}
