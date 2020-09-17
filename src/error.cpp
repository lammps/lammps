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

#include "error.h"

#include "accelerator_kokkos.h"
#include "input.h"
#include "output.h"
#include "universe.h"

#if defined(LAMMPS_EXCEPTIONS)
#include "update.h"
#endif

using namespace LAMMPS_NS;

// helper function to truncate a string to a segment starting with "src/";

static std::string truncpath(const std::string &path)
{
  std::size_t found = path.find("src/");
  if (found != std::string::npos)
    return path.substr(found);
  else return path;
}

/* ---------------------------------------------------------------------- */

Error::Error(LAMMPS *lmp) : Pointers(lmp) {
#ifdef LAMMPS_EXCEPTIONS
  last_error_message.clear();
  last_error_type = ERROR_NONE;
#endif
}

/* ----------------------------------------------------------------------
   called by all procs in universe
   close all output, screen, and log files in world and universe
   no abort, so insure all procs in universe call, else will hang
------------------------------------------------------------------------- */

void Error::universe_all(const std::string &file, int line, const std::string &str)
{
  MPI_Barrier(universe->uworld);
  std::string mesg = fmt::format("ERROR: {} ({}:{})\n",
                                 str,truncpath(file),line);
  if (universe->me == 0) {
    if (universe->uscreen)  fputs(mesg.c_str(),universe->uscreen);
    if (universe->ulogfile) fputs(mesg.c_str(),universe->ulogfile);
  }

  if (output) delete output;
  if (universe->nworlds > 1) {
    if (screen && screen != stdout) fclose(screen);
    if (logfile) fclose(logfile);
  }
  if (universe->ulogfile) fclose(universe->ulogfile);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  // update may be a null pointer when catching command line errors

  if (update) update->whichflag = 0;

  throw LAMMPSException(mesg);
#else
  if (lmp->kokkos) Kokkos::finalize();
  MPI_Finalize();
  exit(1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   forces abort of entire universe if any proc in universe calls
------------------------------------------------------------------------- */

void Error::universe_one(const std::string &file, int line, const std::string &str)
{
  std::string mesg = fmt::format("ERROR on proc {}: {} ({}:{})\n",
                                 universe->me,str,truncpath(file),line);
  if (universe->uscreen) fputs(mesg.c_str(),universe->uscreen);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  // update may be a null pointer when catching command line errors

  if (update) update->whichflag = 0;

  throw LAMMPSAbortException(mesg, universe->uworld);
#else
  MPI_Abort(universe->uworld,1);
  exit(1); // to trick "smart" compilers into believing this does not return
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in universe
   prints a warning message to the screen
------------------------------------------------------------------------- */

void Error::universe_warn(const std::string &file, int line, const std::string &str)
{
  if (universe->uscreen)
    fmt::print(universe->uscreen,"WARNING on proc {}: {} ({}:{})\n",
               universe->me,str,truncpath(file),line);
}

/* ----------------------------------------------------------------------
   called by all procs in one world
   close all output, screen, and log files in world
   insure all procs in world call, else will hang
   force MPI_Abort if running in multi-partition mode
------------------------------------------------------------------------- */

void Error::all(const std::string &file, int line, const std::string &str)
{
  MPI_Barrier(world);

  int me;
  std::string lastcmd = "(unknown)";

  MPI_Comm_rank(world,&me);

  if (me == 0) {
    if (input && input->line) lastcmd = input->line;
    utils::logmesg(lmp,fmt::format("ERROR: {} ({}:{})\nLast command: {}\n",
                                   str,truncpath(file),line,lastcmd));
  }

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  // update may be a null pointer when catching command line errors

  if (update) update->whichflag = 0;

  std::string msg = fmt::format("ERROR: {} ({}:{})\n",
                                str, truncpath(file), line);

  if (universe->nworlds > 1) {
    throw LAMMPSAbortException(msg, universe->uworld);
  }

  throw LAMMPSException(msg);
#else
  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  if (universe->nworlds > 1) MPI_Abort(universe->uworld,1);
  if (lmp->kokkos) Kokkos::finalize();
  MPI_Finalize();
  exit(1);
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in world
   write to world screen only if non-nullptr on this proc
   always write to universe screen
   forces abort of entire world (and universe) if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const std::string &file, int line, const std::string &str)
{
  int me;
  std::string lastcmd = "(unknown)";
  MPI_Comm_rank(world,&me);

  if (input && input->line) lastcmd = input->line;
  std::string mesg = fmt::format("ERROR on proc {}: {} ({}:{})\n",
                                 me,str,truncpath(file),line,lastcmd);
  utils::logmesg(lmp,mesg);

  if (universe->nworlds > 1)
    if (universe->uscreen)
      fputs(mesg.c_str(),universe->uscreen);

#ifdef LAMMPS_EXCEPTIONS

  // allow commands if an exception was caught in a run
  // update may be a null pointer when catching command line errors

  if (update) update->whichflag = 0;

  throw LAMMPSAbortException(mesg, world);
#else
  if (screen) fflush(screen);
  if (logfile) fflush(logfile);
  MPI_Abort(world,1);
  exit(1); // to trick "smart" compilers into believing this does not return
#endif
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-nullptr on this proc since could be file
------------------------------------------------------------------------- */

void Error::warning(const std::string &file, int line, const std::string &str, int logflag)
{
  std::string mesg = fmt::format("WARNING: {} ({}:{})\n",
                                 str,truncpath(file),line);
  if (screen) fputs(mesg.c_str(),screen);
  if (logflag && logfile) fputs(mesg.c_str(),logfile);
}

/* ----------------------------------------------------------------------
   called by one proc in world, typically proc 0
   write message to screen and logfile (if logflag is set)
------------------------------------------------------------------------- */

void Error::message(const std::string &file, int line, const std::string &str, int logflag)
{
  std::string mesg = fmt::format("{} ({}:{})\n",str,truncpath(file),line);

  if (screen) fputs(mesg.c_str(),screen);
  if (logflag && logfile) fputs(mesg.c_str(),logfile);
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

  if (lmp->kokkos) Kokkos::finalize();
  MPI_Finalize();
  exit(status);
}

#ifdef LAMMPS_EXCEPTIONS
/* ----------------------------------------------------------------------
   return the last error message reported by LAMMPS (only used if
   compiled with -DLAMMPS_EXCEPTIONS)
------------------------------------------------------------------------- */

std::string Error::get_last_error() const
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

void Error::set_last_error(const std::string &msg, ErrorType type)
{
  last_error_message = msg;
  last_error_type = type;
}
#endif
