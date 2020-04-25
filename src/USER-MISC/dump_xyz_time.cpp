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

#include "dump_xyz_time.h"
#include <cstring>
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "force.h"

using namespace LAMMPS_NS;

#define ONELINE 128
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpXYZ_Time::DumpXYZ_Time(LAMMPS *lmp, int narg, char **arg) : 
    DumpXYZ(lmp, narg-1, RemoveTimeIncrement(narg, arg))
{
  if (narg != 6) error->all(FLERR,"Wrong number of parameters for dump xyz/time command");
  time_every = force->numeric(FLERR,arg[5]);
  next_time  = 0;
  tol = 1e-3;
  write_time = false;
}

/* ---------------------------------------------------------------------- */

DumpXYZ_Time::~DumpXYZ_Time()
{
}

/* ---------------------------------------------------------------------- */

void DumpXYZ_Time::write_header(bigint n)
{
  update->update_time();
  if (next_time-(update->atime) < time_every*tol) {
      write_time = true;
      next_time +=time_every;
  }
  else write_time = false;
  if (me == 0 && write_time) {
    fprintf(fp,BIGINT_FORMAT "\n",n);
    fprintf(fp,"Atoms. Timestep: " BIGINT_FORMAT " Time: %f" "\n",update->ntimestep, update->atime);
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZ_Time::write_data(int n, double *mybuf)
{
  if (write_time)
    DumpXYZ::write_data(n, mybuf);
}

