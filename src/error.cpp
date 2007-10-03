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

#include "mpi.h"
#include "stdlib.h"
#include "error.h"
#include "universe.h"
#include "memory.h"
#include "output.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Error::Error(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   called by all procs in universe
   close all output, screen, and log files in world and universe
------------------------------------------------------------------------- */

void Error::universe_all(const char *str)
{
  MPI_Barrier(universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,"ERROR: %s\n",str);
    if (universe->ulogfile) fprintf(universe->ulogfile,"ERROR: %s\n",str);
  }

  if (output) delete output;
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in universe
------------------------------------------------------------------------- */

void Error::universe_one(const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"ERROR on proc %d: %s\n",universe->me,str);
  MPI_Abort(universe->uworld,1);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
------------------------------------------------------------------------- */

void Error::all(const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if (screen) fprintf(screen,"ERROR: %s\n",str);
    if (logfile) fprintf(logfile,"ERROR: %s\n",str);
  }

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-NULL on this proc
   always write to universe screen 
------------------------------------------------------------------------- */

void Error::one(const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);
  if (screen) fprintf(screen,"ERROR on proc %d: %s\n",me,str);
  if (universe->nworlds > 1)
    fprintf(universe->uscreen,"ERROR on proc %d: %s\n",universe->me,str);
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file 
------------------------------------------------------------------------- */

void Error::warning(const char *str)
{
  if (screen) fprintf(screen,"WARNING: %s\n",str);
}
