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

#include "dump_custom_time.h"
#include "dump_xyz_time.h" //helper function
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_store.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

#define INVOKED_PERATOM 8
#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCustom_Time::DumpCustom_Time(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg-1, RemoveTimeIncrement(narg, arg))
{
  if (narg == 6) error->all(FLERR,"No dump custom arguments specified");
  time_every = force->numeric(FLERR,arg[5]);
  next_time  = 0;
  tol = 1e-3;
  write_time = 10;
}

/* ---------------------------------------------------------------------- */

DumpCustom_Time::~DumpCustom_Time()
{
}

/* ---------------------------------------------------------------------- */

void DumpCustom_Time::write_header(bigint ndump)
{
  if (write_time == 0)
    DumpCustom::write_header(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustom_Time::write_data(int n, double *mybuf)
{
  if (write_time == 0) 
    (this->*write_choice)(n,mybuf);
}

int DumpCustom_Time::count()
{
  if (write_time >= 0) { //if we are far from write
    update->update_time();
    if (update->atime == 0.0) {
      write_time = 0;
      next_time += time_every;
      return DumpCustom::count();
    }
    else if (next_time-(update->atime+update->dt) < time_every*tol) { //check if it is in next step
      write_time=-1;
      next_time += time_every;
      if (ncompute) {
        for (int i = 0; i < ncompute; i++)
          compute[i]->addstep(update->ntimestep+1);
        }
      return 0;
    }
    else {
      write_time = 1;
      return 0;
    }
  }
  else if (write_time == -1) {
    write_time++;
    return DumpCustom::count();
  }

}
