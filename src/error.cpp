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

#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "error.h"
#include "universe.h"
#include "update.h"
#include "output.h"
#include "input.h"

using namespace LAMMPS_NS;

// helper function to truncate a string to a segment starting with "src/";

static const char *truncpath(const char *path)
{
   if (path) {
     int len = strlen(path);
     for (int i = len-4; i > 0; --i) {
        if (strncmp("src/",path+i,4) == 0)
          return path+i;
     }
   }
   return path;
}

/* ---------------------------------------------------------------------- */

Error::Error(LAMMPS *lmp) : Pointers(lmp) {
#ifdef LAMMPS_EXCEPTIONS
  last_error_message = NULL;
  last_error_type = ERROR_NONE;
#endif
}

/* ----------------------------------------------------------------------
   called by all procs in universe
   close all output, screen, and log files in world and universe
   no abort, so insure all procs in universe call, else will hang
------------------------------------------------------------------------- */

void Error::universe_all(const char *file, int line, const char *str)
{
  MPI_Barrier(universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) fprintf(universe->uscreen,
                                   "ERROR: %s (%s:%d)\n",str,truncpath(file),line);
    if (universe->ulogfile) fprintf(universe->ulogfile,
                                    "ERROR: %s (%s:%d)\n",str,truncpath(file),line);
  }

  if (output) delete output;
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  update->whichflag = 0;

  char msg[100];
  snprintf(msg, 100, "ERROR: %s (%s:%d)\n", str, truncpath(file), line);
  throw LAMMPSException(msg);
#else
  MPI_Finalize();
  exit(1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   forces abort of entire universe if any proc in universe calls
------------------------------------------------------------------------- */

void Error::universe_one(const char *file, int line, const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
            universe->me,str,truncpath(file),line);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  update->whichflag = 0;

  char msg[100];
  snprintf(msg, 100, "ERROR: %s (%s:%d)\n", str, truncpath(file), line);
  throw LAMMPSAbortException(msg, universe->uworld);
#else
  MPI_Abort(universe->uworld,1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   prints a warning message to the screen
------------------------------------------------------------------------- */

void Error::universe_warn(const char *file, int line, const char *str)
{
  if (universe->uscreen)
    fprintf(universe->uscreen,"WARNING on proc %d: %s (%s:%d)\n",
            universe->me,str,truncpath(file),line);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   insure all procs in world call, else will hang
   force MPI_Abort if running in multi-partition mode
------------------------------------------------------------------------- */

void Error::all(const char *file, int line, const char *str)
{
  MPI_Barrier(world);

  int me;
  const char *lastcmd = (const char*)"(unknown)";

  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if (input && input->line) lastcmd = input->line;
    if (screen) fprintf(screen,"ERROR: %s (%s:%d)\n"
                        "Last command: %s\n",
                        str,truncpath(file),line,lastcmd);
    if (logfile) fprintf(logfile,"ERROR: %s (%s:%d)\n"
                         "Last command: %s\n",
                         str,truncpath(file),line,lastcmd);
  }

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  update->whichflag = 0;

  char msg[100];
  snprintf(msg, 100, "ERROR: %s (%s:%d)\n", str, truncpath(file), line);

  if (universe->nworlds > 1) {
    throw LAMMPSAbortException(msg, universe->uworld);
  }

  throw LAMMPSException(msg);
#else
  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  MPI_Finalize();
  exit(1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-NULL on this proc
   always write to universe screen
   forces abort of entire world (and universe) if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const char *file, int line, const char *str)
{
  int me;
  const char *lastcmd = (const char*)"(unknown)";
  MPI_Comm_rank(world,&me);

  if (input && input->line) lastcmd = input->line;
  if (screen) fprintf(screen,"ERROR on proc %d: %s (%s:%d)\n"
                      "Last command: %s\n",
                      me,str,truncpath(file),line,lastcmd);
  if (logfile) fprintf(logfile,"ERROR on proc %d: %s (%s:%d)\n"
                       "Last command: %s\n",
                       me,str,truncpath(file),line,lastcmd);

  if (universe->nworlds > 1)
    if (universe->uscreen)
      fprintf(universe->uscreen,"ERROR on proc %d: %s (%s:%d)\n",
              universe->me,str,truncpath(file),line);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  update->whichflag = 0;

  char msg[100];
  snprintf(msg, 100, "ERROR on proc %d: %s (%s:%d)\n", me, str, truncpath(file), line);
  throw LAMMPSAbortException(msg, world);
#else
  MPI_Abort(world,1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file
------------------------------------------------------------------------- */

void Error::warning(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"WARNING: %s (%s:%d)\n",str,truncpath(file),line);
  if (logflag && logfile) fprintf(logfile,"WARNING: %s (%s:%d)\n",
                                  str,truncpath(file),line);
}

/* ----------------------------------------------------------------------
   called by one proc in world, typically proc 0
   write message to screen and logfile (if logflag is set)
------------------------------------------------------------------------- */

void Error::message(const char *file, int line, const char *str, int logflag)
{
  if (screen) fprintf(screen,"%s (%s:%d)\n",str,truncpath(file),line);
  if (logflag && logfile) fprintf(logfile,"%s (%s:%d)\n",str,truncpath(file),line);
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   called by all procs in one world
   close all output, screen, and log files in world
   no abort, so insure all procs in world call, else will hang
------------------------------------------------------------------------- */

void Error::done(int status)
{
  MPI_Barrier(world);

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  MPI_Finalize();
  exit(status);
}

#ifdef LAMMPS_EXCEPTIONS
/* ----------------------------------------------------------------------
   return the last error message reported by LAMMPS (only used if
   compiled with -DLAMMPS_EXCEPTIONS)
------------------------------------------------------------------------- */

char * Error::get_last_error() const
{
  return last_error_message;
}

/* ----------------------------------------------------------------------
   return the type of the last error reported by LAMMPS (only used if
   compiled with -DLAMMPS_EXCEPTIONS)
------------------------------------------------------------------------- */

ErrorType Error::get_last_error_type() const
{
  return last_error_type;
}

/* ----------------------------------------------------------------------
   set the last error message and error type
   (only used if compiled with -DLAMMPS_EXCEPTIONS)
------------------------------------------------------------------------- */

void Error::set_last_error(const char * msg, ErrorType type)
{
  delete [] last_error_message;

  if(msg) {
    last_error_message = new char[strlen(msg)+1];
    strcpy(last_error_message, msg);
  } else {
    last_error_message = NULL;
  }
  last_error_type = type;
}
#endif
