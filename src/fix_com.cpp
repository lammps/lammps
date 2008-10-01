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
#include "fix_com.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixCOM::FixCOM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal fix com command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix com command");
  first = 1;

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix com file %s",arg[4]);
      error->one(str);
    }
  }

  if (me == 0) {
    fprintf(fp,"# Center-of-mass for group %s\n",group->names[igroup]);
    fprintf(fp,"# TimeStep x y z\n");
  }
}

/* ---------------------------------------------------------------------- */

FixCOM::~FixCOM()
{
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixCOM::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCOM::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixCOM::setup(int vflag)
{
  if (first) end_of_step();
  first = 0;
}

/* ---------------------------------------------------------------------- */

void FixCOM::end_of_step()
{
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  if (me == 0) {
    fprintf(fp,"%d %g %g %g\n",update->ntimestep,xcm[0],xcm[1],xcm[2]);
    fflush(fp);
  }
}
