// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Michal Kanski (Jagiellonian U) for simulation time dumps
------------------------------------------------------------------------- */

#include "output.h"
#include "style_dump.h"         // IWYU pragma: keep

#include "comm.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "thermo.h"
#include "update.h"
#include "variable.h"
#include "write_restart.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define DELTA 1
#define EPSDT 1.0e-6

enum {SETUP, WRITE, RESET_DT};

/* ----------------------------------------------------------------------
   one instance per dump style in style_dump.h
------------------------------------------------------------------------- */

template <typename T> static Dump *dump_creator(LAMMPS *lmp, int narg, char ** arg)
{
  return new T(lmp, narg, arg);
}

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(LAMMPS *lmp) : Pointers(lmp)
{
  // create default computes for temp,pressure,pe

  modify->add_compute("thermo_temp all temp");
  modify->add_compute("thermo_press all pressure thermo_temp");
  modify->add_compute("thermo_pe all pe");

  // create default Thermo class

  auto newarg = new char*[1];
  newarg[0] = (char *) "one";
  thermo = new Thermo(lmp,1,newarg);
  delete[] newarg;

  thermo_every = 0;
  var_thermo = nullptr;

  ndump = 0;
  max_dump = 0;
  any_time_dumps = 0;
  next_dump_any = next_time_dump_any = MAXBIGINT;
  mode_dump = nullptr;
  every_dump = nullptr;
  every_time_dump = nullptr;
  next_dump = nullptr;
  next_time_dump = nullptr;
  last_dump = nullptr;
  var_dump = nullptr;
  ivar_dump = nullptr;
  dump = nullptr;

  restart_flag = restart_flag_single = restart_flag_double = 0;
  restart_every_single = restart_every_double = 0;
  last_restart = -1;
  restart1 = restart2a = restart2b = nullptr;
  var_restart_single = var_restart_double = nullptr;
  restart = nullptr;

  dump_map = new DumpCreatorMap();

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  (*dump_map)[#key] = &dump_creator<Class>;
#include "style_dump.h"         // IWYU pragma: keep
#undef DumpStyle
#undef DUMP_CLASS
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Output::~Output()
{
  if (thermo) delete thermo;
  delete[] var_thermo;

  memory->destroy(mode_dump);
  memory->destroy(every_dump);
  memory->destroy(every_time_dump);
  memory->destroy(next_dump);
  memory->destroy(next_time_dump);
  memory->destroy(last_dump);
  for (int i = 0; i < ndump; i++) delete[] var_dump[i];
  memory->sfree(var_dump);
  memory->destroy(ivar_dump);
  for (int i = 0; i < ndump; i++) delete dump[i];
  memory->sfree(dump);

  delete[] restart1;
  delete[] restart2a;
  delete[] restart2b;
  delete[] var_restart_single;
  delete[] var_restart_double;
  delete restart;

  delete dump_map;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  thermo->init();
  if (var_thermo) {
    ivar_thermo = input->variable->find(var_thermo);
    if (ivar_thermo < 0)
      error->all(FLERR,"Variable name for thermo every does not exist");
    if (!input->variable->equalstyle(ivar_thermo))
      error->all(FLERR,"Variable for thermo every is invalid style");
  }

  for (int i = 0; i < ndump; i++) dump[i]->init();
  any_time_dumps = 0;
  for (int i = 0; i < ndump; i++) {
    if (mode_dump[i]) any_time_dumps = 1;
    if ((mode_dump[i] == 0 && every_dump[i] == 0) ||
        (mode_dump[i] == 1 && every_time_dump[i] == 0.0)) {
      ivar_dump[i] = input->variable->find(var_dump[i]);
      if (ivar_dump[i] < 0)
        error->all(FLERR,"Variable name for dump every or delta does not exist");
      if (!input->variable->equalstyle(ivar_dump[i]))
        error->all(FLERR,"Variable for dump every or delta is invalid style");
    }
  }

  if (restart_flag_single && restart_every_single == 0) {
    ivar_restart_single = input->variable->find(var_restart_single);
    if (ivar_restart_single < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equalstyle(ivar_restart_single))
      error->all(FLERR,"Variable for restart is invalid style");
  }
  if (restart_flag_double && restart_every_double == 0) {
    ivar_restart_double = input->variable->find(var_restart_double);
    if (ivar_restart_double < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equalstyle(ivar_restart_double))
      error->all(FLERR,"Variable for restart is invalid style");
  }
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do thermo last, so will print after memory_usage
   memflag = 0/1 for printing out memory usage
------------------------------------------------------------------------- */

void Output::setup(int memflag)
{
  bigint ntimestep = update->ntimestep;

  // print memory usage unless being called between multiple runs

  if (memflag) memory_usage();

  // set next_thermo to multiple of every or variable eval if var defined
  // ensure thermo output on last step of run
  // thermo may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  thermo->header();
  thermo->compute(0);
  last_thermo = ntimestep;

  if (var_thermo) {
    next_thermo = static_cast<bigint>
      (input->variable->compute_equal(ivar_thermo));
    if (next_thermo <= ntimestep)
      error->all(FLERR,"Thermo every variable returned a bad timestep");
  } else if (thermo_every) {
    next_thermo = (ntimestep/thermo_every)*thermo_every + thermo_every;
    next_thermo = MIN(next_thermo,update->laststep);
  } else next_thermo = update->laststep;

  modify->addstep_compute(next_thermo);

  // consider all dumps
  // decide whether to write snapshot and/or calculate next step for dump

  if (ndump && update->restrict_output == 0) {
    next_dump_any = next_time_dump_any = MAXBIGINT;

    for (int idump = 0; idump < ndump; idump++) {

      // wrap step dumps that invoke computes or do variable eval with clear/add
      // see NOTE in write() about also wrapping time dumps

      if (mode_dump[idump] == 0 && (dump[idump]->clearstep || var_dump[idump]))
        modify->clearstep_compute();

      // write a snapshot at setup only if any of these 3 conditions hold
      // (1) this is first run since dump was created and its first_flag = 0
      // (2) mode_dump = 0 and timestep is multiple of every_dump
      // (3) mode_dump = 1 and time is multiple of every_time_dump (within EPSDT)
      // (2) and (3) only apply for non-variable dump intervals
      // finally, do not write if same snapshot written previously,
      //   i.e. on last timestep of previous run

      int writeflag = 0;

      if (last_dump[idump] < 0 && dump[idump]->first_flag == 1) writeflag = 1;

      if (mode_dump[idump] == 0) {
        if (every_dump[idump] && (ntimestep % every_dump[idump] == 0))
          writeflag = 1;
      } else {
        if (every_time_dump[idump] > 0.0) {
          double tcurrent = update->atime +
            (ntimestep - update->atimestep) * update->dt;
          double remainder = fmod(tcurrent,every_time_dump[idump]);
          if ((remainder < EPSDT*update->dt) ||
              (every_time_dump[idump] - remainder < EPSDT*update->dt))
            writeflag = 1;
        }
      }

      if (last_dump[idump] == ntimestep) writeflag = 0;

      // perform dump

      if (writeflag) {
        dump[idump]->write();
        last_dump[idump] = ntimestep;
      }

      // calculate timestep or time for next dump
      // set next_dump and next_time_dump

      calculate_next_dump(SETUP,idump,ntimestep);

      // if dump not written now, use addstep_compute_all()
      // since don't know what computes the dump will invoke

      if (mode_dump[idump] == 0 && (dump[idump]->clearstep || var_dump[idump])) {
        if (writeflag) modify->addstep_compute(next_dump[idump]);
        else modify->addstep_compute_all(next_dump[idump]);
      }

      if (mode_dump[idump] && (dump[idump]->clearstep || var_dump[idump]))
        next_time_dump_any = MIN(next_time_dump_any,next_dump[idump]);
      next_dump_any = MIN(next_dump_any,next_dump[idump]);
    }

    // if no dumps, set next_dump_any to last+1 so will not influence next

  } else next_dump_any = update->laststep + 1;

  // do not write restart files at start of run
  // set next_restart values to multiple of every or variable value
  // wrap variable eval with clear/add
  // if no restarts, set next_restart to last+1 so will not influence next

  if (restart_flag && update->restrict_output == 0) {
    if (restart_flag_single) {
      if (restart_every_single) {
        next_restart_single =
          (ntimestep/restart_every_single)*restart_every_single +
          restart_every_single;
      } else {
        auto  nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad next timestep: {} vs {}",
                     nextrestart, ntimestep);
        next_restart_single = nextrestart;
      }
    } else next_restart_single = update->laststep + 1;
    if (restart_flag_double) {
      if (restart_every_double)
        next_restart_double =
          (ntimestep/restart_every_double)*restart_every_double +
          restart_every_double;
      else {
        auto  nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad next timestep: {} vs {}",
                     nextrestart, ntimestep);
        next_restart_double = nextrestart;
      }
    } else next_restart_double = update->laststep + 1;
    next_restart = MIN(next_restart_single,next_restart_double);
  } else next_restart = update->laststep + 1;

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last output doesn't
   do dump/restart before thermo so thermo CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  // ensure next_thermo forces output on last step of run
  // thermo may invoke computes so wrap with clear/add

  if (next_thermo == ntimestep) {
    modify->clearstep_compute();
    if (last_thermo != ntimestep) thermo->compute(1);
    last_thermo = ntimestep;
    if (var_thermo) {
      next_thermo = static_cast<bigint>
        (input->variable->compute_equal(ivar_thermo));
      if (next_thermo <= ntimestep)
        error->all(FLERR,"Thermo every variable returned a bad timestep");
    } else if (thermo_every) next_thermo += thermo_every;
    else next_thermo = update->laststep;
    next_thermo = MIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  }

  // perform dump if its next_dump = current ntimestep
  //   but not if it was already written on this step
  // set next_dump and also next_time_dump for mode_dump = 1
  // set next_dump_any to smallest next_dump
  // wrap step dumps that invoke computes or do variable eval with clear/add
  // NOTE:
  //   not wrapping time dumps means that Integrate::ev_set()
  //     needs to trigger all per-atom eng/virial computes
  //     on a timestep where any time dump will be output
  //   could wrap time dumps as well, if timestep size did not vary
  //   if wrap when timestep size varies frequently,
  //     then can do many unneeded addstep() --> inefficient
  //   hard to know if timestep varies, since run every could change it
  //   can't remove an uneeded addstep from a compute, b/c don't know
  //     what other command may have added it

  if (next_dump_any == ntimestep) {
    next_dump_any = next_time_dump_any = MAXBIGINT;

    for (int idump = 0; idump < ndump; idump++) {

      if (next_dump[idump] == ntimestep) {
        if (last_dump[idump] == ntimestep) continue;

        if (mode_dump[idump] == 0 &&
            (dump[idump]->clearstep || var_dump[idump]))
          modify->clearstep_compute();

        // perform dump
        // set next_dump and next_time_dump

        dump[idump]->write();
        last_dump[idump] = ntimestep;
        calculate_next_dump(WRITE,idump,ntimestep);

        if (mode_dump[idump] == 0 &&
            (dump[idump]->clearstep || var_dump[idump]))
          modify->addstep_compute(next_dump[idump]);
      }

      if (mode_dump[idump] && (dump[idump]->clearstep || var_dump[idump]))
        next_time_dump_any = MIN(next_time_dump_any,next_dump[idump]);
      next_dump_any = MIN(next_dump_any,next_dump[idump]);
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // next restart variable may invoke computes so wrap with clear/add

  if (next_restart == ntimestep) {
    if (next_restart_single == ntimestep) {
      std::string file = restart1;
      std::size_t found = file.find('*');
      if (found != std::string::npos)
        file.replace(found,1,fmt::format("{}",update->ntimestep));

      if (last_restart != ntimestep) restart->write(file);

      if (restart_every_single) next_restart_single += restart_every_single;
      else {
        modify->clearstep_compute();
        auto  nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad next timestep: {} vs {}",
                     nextrestart, ntimestep);
        next_restart_single = nextrestart;
        modify->addstep_compute(next_restart_single);
      }
    }

    if (next_restart_double == ntimestep) {
      if (last_restart != ntimestep) {
        if (restart_toggle == 0) {
          restart->write(restart2a);
          restart_toggle = 1;
        } else {
          restart->write(restart2b);
          restart_toggle = 0;
        }
      }

      if (restart_every_double) next_restart_double += restart_every_double;
      else {
        modify->clearstep_compute();
        auto  nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad next timestep: {} <= {}",
                     nextrestart, ntimestep);
        next_restart_double = nextrestart;
        modify->addstep_compute(next_restart_double);
      }
    }
    last_restart = ntimestep;
    next_restart = MIN(next_restart_single,next_restart_double);
  }

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   force a snapshot to be written for all dumps
   called from PRD and TAD
------------------------------------------------------------------------- */

void Output::write_dump(bigint ntimestep)
{
  for (int idump = 0; idump < ndump; idump++) {
    dump[idump]->write();
    last_dump[idump] = ntimestep;
  }
}

/* ----------------------------------------------------------------------
   calculate when next dump occurs for Dump instance idump
   operates in one of two modes, based on mode_dump flag
   for timestep mode, set next_dump
   for simulation time mode, set next_time_dump and next_dump
   which flag depends on caller
   SETUP = from setup() at start of run
   WRITE = from write() during run each time a dump file is written
   RESET_DT = from reset_dt() called from fix dt/reset when it changes timestep size
------------------------------------------------------------------------- */

void Output::calculate_next_dump(int which, int idump, bigint ntimestep)
{
  // dump mode is by timestep
  // just set next_dump

  if (mode_dump[idump] == 0) {

    if (every_dump[idump]) {

      // which = SETUP: next_dump = next multiple of every_dump
      // which = WRITE: increment next_dump by every_dump
      //                current step is already multiple of every_dump

      if (which == SETUP)
        next_dump[idump] = (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      else if (which == WRITE)
        next_dump[idump] += every_dump[idump];

    } else {
      next_dump[idump] = static_cast<bigint>(input->variable->compute_equal(ivar_dump[idump]));
      if (next_dump[idump] <= ntimestep)
        error->all(FLERR,"Dump {} every variable {} returned a bad timestep: {}",
                   dump[idump]->id, var_dump[idump], next_dump[idump]);
    }

    // dump mode is by simulation time
    // set next_time_dump and next_dump

  } else {

    bigint nextdump;
    double nexttime;
    double tcurrent = update->atime + (ntimestep - update->atimestep) * update->dt;

    if (every_time_dump[idump] > 0.0) {

      // which = SETUP: nexttime = next multiple of every_time_dump
      // which = WRITE: increment nexttime by every_time_dump
      // which = RESET_DT: no change to previous nexttime (only timestep has changed)

      switch (which) {
      case SETUP:
        nexttime = static_cast<bigint> (tcurrent/every_time_dump[idump]) *
          every_time_dump[idump] + every_time_dump[idump];
        break;

      case WRITE:
        nexttime = next_time_dump[idump] + every_time_dump[idump];
        break;

      case RESET_DT:
        nexttime = next_time_dump[idump];
        break;

      default:
        nexttime = 0;
        error->all(FLERR,"Unexpected argument to calculate_next_dump");
      }

      nextdump = ntimestep +
        static_cast<bigint> ((nexttime - tcurrent - EPSDT*update->dt) / update->dt) + 1;

      // if delta is too small to reach next timestep, use multiple of delta

      if (nextdump == ntimestep) {
        double tnext = update->atime + (ntimestep + 1 - update->atimestep) * update->dt;
        int multiple = static_cast<int>((tnext - nexttime) / every_time_dump[idump]);
        nexttime = nexttime + (multiple + 1) * every_time_dump[idump];
        nextdump = ntimestep +
          static_cast<bigint> ((nexttime - tcurrent - EPSDT*update->dt) / update->dt) + 1;
      }

    } else {

      // do not re-evaulate variable for which = RESET_DT, leave nexttime as-is
      // unless next_time_dump < 0.0, which means variable never yet evaluated

      if (which < RESET_DT || next_time_dump[idump] < 0.0) {
        nexttime = input->variable->compute_equal(ivar_dump[idump]);
      } else
        nexttime = next_time_dump[idump];

      if (nexttime <= tcurrent)
        error->all(FLERR,"Dump every/time variable returned a bad time");

      nextdump = ntimestep +
        static_cast<bigint> ((nexttime - tcurrent - EPSDT*update->dt) / update->dt) + 1;
      if (nextdump <= ntimestep)
        error->all(FLERR,"Dump every/time variable too small for next timestep");
    }

    next_time_dump[idump] = nexttime;
    next_dump[idump] = nextdump;
  }
}

/* ---------------------------------------------------------------------- */

int Output::check_time_dumps(bigint ntimestep)
{
  int nowflag = 0;
  for (int i = 0; i < ndump; i++)
    if (mode_dump[i] && next_dump[i] == ntimestep) nowflag = 1;

  return nowflag;
}

/* ----------------------------------------------------------------------
   force restart file(s) to be written
   called from PRD and TAD
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_flag_single) {
    std::string file = restart1;
    std::size_t found = file.find('*');
    if (found != std::string::npos)
      file.replace(found,1,fmt::format("{}",update->ntimestep));
    restart->write(file);
  }

  if (restart_flag_double) {
    if (restart_toggle == 0) {
      restart->write(restart2a);
      restart_toggle = 1;
    } else {
      restart->write(restart2b);
      restart_toggle = 0;
    }
  }

  last_restart = ntimestep;
}

/* ----------------------------------------------------------------------
   timestep is being changed, called by update->reset_timestep()
   for dumps, require that no dump is "active"
   meaning that a snapshot has already been output
   reset next output values for restart and thermo
   reset to smallest value >= new timestep
   if next timestep set by variable evaluation,
   eval for ntimestep-1, so current ntimestep can be returned if needed
   no guarantee that variable can be evaluated for ntimestep-1
   e.g. if it depends on computes, but live with that rare case for now
------------------------------------------------------------------------- */

void Output::reset_timestep(bigint ntimestep)
{
  next_dump_any = MAXBIGINT;
  for (int idump = 0; idump < ndump; idump++)
    if ((last_dump[idump] >= 0) && !update->whichflag && !dump[idump]->multifile)
      error->all(FLERR, "Cannot reset timestep with active dump - must undump first");

  if (restart_flag_single) {
    if (restart_every_single) {
      next_restart_single =
        (ntimestep/restart_every_single)*restart_every_single;
      if (next_restart_single < ntimestep)
        next_restart_single += restart_every_single;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      auto  nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_single));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad next timestep: {} <= {}",
                   nextrestart, ntimestep);
      update->ntimestep++;
      next_restart_single = nextrestart;
      modify->addstep_compute(next_restart_single);
    }
  } else next_restart_single = update->laststep + 1;

  if (restart_flag_double) {
    if (restart_every_double) {
      next_restart_double =
        (ntimestep/restart_every_double)*restart_every_double;
      if (next_restart_double < ntimestep)
        next_restart_double += restart_every_double;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      auto  nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_double));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad next timestep: {} <= {}",
                   nextrestart, ntimestep);
      update->ntimestep++;
      next_restart_double = nextrestart;
      modify->addstep_compute(next_restart_double);
    }
  } else next_restart_double = update->laststep + 1;

  next_restart = MIN(next_restart_single,next_restart_double);

  if (var_thermo) {
    modify->clearstep_compute();
    update->ntimestep--;
    next_thermo = static_cast<bigint>
      (input->variable->compute_equal(ivar_thermo));
    if (next_thermo < ntimestep)
      error->all(FLERR,"Thermo_modify every variable returned a bad timestep");
    update->ntimestep++;
    next_thermo = MIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  } else if (thermo_every) {
    next_thermo = (ntimestep/thermo_every)*thermo_every;
    if (next_thermo < ntimestep) next_thermo += thermo_every;
    next_thermo = MIN(next_thermo,update->laststep);
  } else next_thermo = update->laststep;

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   timestep size is being changed
   reset next output values for dumps which have mode_dump=1
   called by fix dt/reset (at end of step)
   or called by timestep command via run every (also at end of step)
------------------------------------------------------------------------- */

void Output::reset_dt()
{
  bigint ntimestep = update->ntimestep;

  next_time_dump_any = MAXBIGINT;

  for (int idump = 0; idump < ndump; idump++) {
    if (mode_dump[idump] == 0) continue;

    // reset next_dump but do not change next_time_dump, 2 arg for reset_dt()
    // do not invoke for a dump already scheduled for this step
    //   since timestep change affects next step

    if (next_dump[idump] != ntimestep)
      calculate_next_dump(RESET_DT,idump,update->ntimestep);

    if (dump[idump]->clearstep || var_dump[idump])
      next_time_dump_any = MIN(next_time_dump_any,next_dump[idump]);
  }

  next_dump_any = MIN(next_dump_any,next_time_dump_any);
  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

Dump *Output::add_dump(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) error->all(FLERR,"Reuse of dump ID: {}", arg[0]);

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find dump group ID: {}", arg[1]);
  if (utils::inumeric(FLERR,arg[3],false,lmp) <= 0)
    error->all(FLERR,"Invalid dump frequency {}", arg[3]);

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **) memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    memory->grow(mode_dump,max_dump,"output:mode_dump");
    memory->grow(every_dump,max_dump,"output:every_dump");
    memory->grow(every_time_dump,max_dump,"output:every_time_dump");
    memory->grow(next_dump,max_dump,"output:next_dump");
    memory->grow(next_time_dump,max_dump,"output:next_time_dump");
    memory->grow(last_dump,max_dump,"output:last_dump");
    var_dump = (char **) memory->srealloc(var_dump,max_dump*sizeof(char *),"output:var_dump");
    memory->grow(ivar_dump,max_dump,"output:ivar_dump");
  }

  // create the Dump
  int idump = ndump;

  if (dump_map->find(arg[2]) != dump_map->end()) {
    DumpCreator &dump_creator = (*dump_map)[arg[2]];
    dump[idump] = dump_creator(lmp, narg, arg);
  } else error->all(FLERR,utils::check_packages_for_style("dump",arg[2],lmp));

  // initialize per-dump data to suitable default values

  mode_dump[idump] = 0;
  every_dump[idump] = utils::inumeric(FLERR,arg[3],false,lmp);
  if (every_dump[idump] <= 0) error->all(FLERR,"Illegal dump command");
  every_time_dump[idump] = 0.0;
  next_time_dump[idump] = -1.0;
  last_dump[idump] = -1;
  var_dump[idump] = nullptr;
  ivar_dump[idump] = -1;
  next_dump[idump] = 0;

  ndump++;
  dump_list = std::vector<Dump *>(dump, dump + ndump);
  return dump[idump];
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify",error);

  // find which dump it is

  auto idump = get_dump_by_id(arg[0]);
  if (!idump) error->all(FLERR,"Could not find dump_modify ID: {}", arg[0]);
  idump->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Dump from list of Dumps
------------------------------------------------------------------------- */

void Output::delete_dump(const std::string &id)
{
  // find which dump it is and delete it

  int idump;
  for (idump = 0; idump < ndump; idump++) if (id == dump[idump]->id) break;
  if (idump == ndump) error->all(FLERR,"Could not find undump ID: {}", id);

  delete dump[idump];
  delete[] var_dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    mode_dump[i-1] = mode_dump[i];
    every_dump[i-1] = every_dump[i];
    every_time_dump[i-1] = every_time_dump[i];
    next_dump[i-1] = next_dump[i];
    next_time_dump[i-1] = next_time_dump[i];
    last_dump[i-1] = last_dump[i];
    var_dump[i-1] = var_dump[i];
    ivar_dump[i-1] = ivar_dump[i];
  }
  ndump--;
  dump[ndump] = nullptr;
  var_dump[ndump] = nullptr;
  dump_list = std::vector<Dump *>(dump, dump + ndump);
}

/* ----------------------------------------------------------------------
   find a dump by ID
   return pointer to dump
------------------------------------------------------------------------- */

Dump *Output::get_dump_by_id(const std::string &id) const
{
  if (id.empty()) return nullptr;
  for (int idump = 0; idump < ndump; idump++) if (id == dump[idump]->id) return dump[idump];
  return nullptr;
}

/* ----------------------------------------------------------------------
   return list of dumps as vector
------------------------------------------------------------------------- */

const std::vector<Dump *> &Output::get_dump_list()
{
  dump_list = std::vector<Dump *>(dump, dump + ndump);
  return dump_list;
}

/* ----------------------------------------------------------------------
   set thermo output frequency from input script
------------------------------------------------------------------------- */

void Output::set_thermo(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal thermo command");

  // always reset var_thermo, so it is possible to switch back from
  // variable spaced thermo outputs to constant spaced ones.

  delete[] var_thermo;
  var_thermo = nullptr;

  if (utils::strmatch(arg[0],"^v_")) {
    var_thermo = utils::strdup(arg[0]+2);
  } else {
    thermo_every = utils::inumeric(FLERR,arg[0],false,lmp);
    if (thermo_every < 0) error->all(FLERR,"Illegal thermo output frequency {}", thermo_every);
  }
}

/* ----------------------------------------------------------------------
   new Thermo style
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "thermo_style", error);

  // don't allow this so that dipole style can safely allocate inertia vector

  if (domain->box_exist == 0)
    error->all(FLERR,"Thermo_style command before simulation box is defined");

  // warn if previous thermo had been modified via thermo_modify command

  if (thermo->modified && comm->me == 0)
    error->warning(FLERR,"New thermo_style command, previous thermo_modify settings will be lost");

  // set thermo = nullptr in case new Thermo throws an error

  delete thermo;
  thermo = nullptr;
  thermo = new Thermo(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   setup restart capability for single or double output files
   if only one filename and it contains no "*", then append ".*"
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "restart", error);

  int every = 0;
  int varflag = 0;

  if (utils::strmatch(arg[0],"^v_")) varflag = 1;
  else every = utils::inumeric(FLERR,arg[0],false,lmp);

  if (!varflag && every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");

    restart_flag = restart_flag_single = restart_flag_double = 0;
    last_restart = -1;

    delete restart;
    restart = nullptr;
    delete[] restart1;
    delete[] restart2a;
    delete[] restart2b;
    restart1 = restart2a = restart2b = nullptr;
    delete[] var_restart_single;
    delete[] var_restart_double;
    var_restart_single = var_restart_double = nullptr;

    return;
  }

  if (narg < 2) error->all(FLERR,"Illegal restart command");

  int nfile = 0;
  if (narg % 2 == 0) nfile = 1;
  else nfile = 2;

  if (nfile == 1) {
    restart_flag = restart_flag_single = 1;

    if (varflag) {
      delete[] var_restart_single;
      var_restart_single = utils::strdup(arg[0]+2);
      restart_every_single = 0;
    } else restart_every_single = every;

    int n = strlen(arg[1]) + 3;
    delete[] restart1;
    restart1 = new char[n];
    strcpy(restart1,arg[1]);
    if (strchr(restart1,'*') == nullptr) strcat(restart1,".*");
  }

  if (nfile == 2) {
    restart_flag = restart_flag_double = 1;

    if (varflag) {
      delete[] var_restart_double;
      var_restart_double = utils::strdup(arg[0]+2);
      restart_every_double = 0;
    } else restart_every_double = every;

    delete[] restart2a;
    delete[] restart2b;
    restart_toggle = 0;
    restart2a = utils::strdup(arg[1]);
    restart2b = utils::strdup(arg[2]);
  }

  // check for multiproc output and an MPI-IO filename
  // if 2 filenames, must be consistent

  int multiproc;
  if (strchr(arg[1],'%')) multiproc = comm->nprocs;
  else multiproc = 0;
  if (nfile == 2) {
    if (multiproc && !strchr(arg[2],'%'))
      error->all(FLERR,"Both restart files must use % or neither");
    if (!multiproc && strchr(arg[2],'%'))
      error->all(FLERR,"Both restart files must use % or neither");
  }

  int mpiioflag;
  if (utils::strmatch(arg[1],"\\.mpiio$")) mpiioflag = 1;
  else mpiioflag = 0;
  if (nfile == 2) {
    if (mpiioflag && !utils::strmatch(arg[2],"\\.mpiio$"))
      error->all(FLERR,"Both restart files must use MPI-IO or neither");
    if (!mpiioflag && utils::strmatch(arg[2],"\\.mpiio$"))
      error->all(FLERR,"Both restart files must use MPI-IO or neither");
  }

  // setup output style and process optional args

  delete restart;
  restart = new WriteRestart(lmp);
  int iarg = nfile+1;
  restart->multiproc_options(multiproc,mpiioflag,narg-iarg,&arg[iarg]);
}

/* ----------------------------------------------------------------------
   sum and print memory usage
   result is only memory on proc 0, not averaged across procs
------------------------------------------------------------------------- */

void Output::memory_usage()
{
  double meminfo[3];
  Info info(lmp);

  info.get_memory_info(meminfo);
  double mbytes = meminfo[0];
  double mbmin,mbavg,mbmax;
  MPI_Reduce(&mbytes,&mbavg,1,MPI_DOUBLE,MPI_SUM,0,world);
  MPI_Reduce(&mbytes,&mbmin,1,MPI_DOUBLE,MPI_MIN,0,world);
  MPI_Reduce(&mbytes,&mbmax,1,MPI_DOUBLE,MPI_MAX,0,world);
  mbavg /= comm->nprocs;

  if (comm->me == 0)
    utils::logmesg(lmp,"Per MPI rank memory allocation (min/avg/max) = "
                   "{:.4} | {:.4} | {:.4} Mbytes\n",mbmin,mbavg,mbmax);
}
