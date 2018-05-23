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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "output.h"
#include "style_dump.h"
#include "atom.h"
#include "neighbor.h"
#include "input.h"
#include "variable.h"
#include "comm.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "thermo.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "dump.h"
#include "write_restart.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 1

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(LAMMPS *lmp) : Pointers(lmp)
{
  // create default computes for temp,pressure,pe

  char **newarg = new char*[4];
  newarg[0] = (char *) "thermo_temp";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);

  newarg[0] = (char *) "thermo_press";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "thermo_temp";
  modify->add_compute(4,newarg);

  newarg[0] = (char *) "thermo_pe";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);

  delete [] newarg;

  // create default Thermo class

  newarg = new char*[1];
  newarg[0] = (char *) "one";
  thermo = new Thermo(lmp,1,newarg);
  delete [] newarg;

  thermo_every = 0;
  var_thermo = NULL;

  ndump = 0;
  max_dump = 0;
  every_dump = NULL;
  next_dump = NULL;
  last_dump = NULL;
  var_dump = NULL;
  ivar_dump = NULL;
  dump = NULL;

  restart_flag = restart_flag_single = restart_flag_double = 0;
  restart_every_single = restart_every_double = 0;
  last_restart = -1;
  restart1 = restart2a = restart2b = NULL;
  var_restart_single = var_restart_double = NULL;
  restart = NULL;

  dump_map = new DumpCreatorMap();

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  (*dump_map)[#key] = &dump_creator<Class>;
#include "style_dump.h"
#undef DumpStyle
#undef DUMP_CLASS
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Output::~Output()
{
  if (thermo) delete thermo;
  delete [] var_thermo;

  memory->destroy(every_dump);
  memory->destroy(next_dump);
  memory->destroy(last_dump);
  for (int i = 0; i < ndump; i++) delete [] var_dump[i];
  memory->sfree(var_dump);
  memory->destroy(ivar_dump);
  for (int i = 0; i < ndump; i++) delete dump[i];
  memory->sfree(dump);

  delete [] restart1;
  delete [] restart2a;
  delete [] restart2b;
  delete [] var_restart_single;
  delete [] var_restart_double;
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
  for (int i = 0; i < ndump; i++)
    if (every_dump[i] == 0) {
      ivar_dump[i] = input->variable->find(var_dump[i]);
      if (ivar_dump[i] < 0)
        error->all(FLERR,"Variable name for dump every does not exist");
      if (!input->variable->equalstyle(ivar_dump[i]))
        error->all(FLERR,"Variable for dump every is invalid style");
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

  // perform dump at start of run only if:
  //   current timestep is multiple of every and last dump not >= this step
  //   this is first run after dump created and firstflag is set
  //   note that variable freq will not write unless triggered by firstflag
  // set next_dump to multiple of every or variable value
  // set next_dump_any to smallest next_dump
  // wrap dumps that invoke computes and variable eval with clear/add
  // if dump not written now, use addstep_compute_all() since don't know
  //   what computes the dump write would invoke
  // if no dumps, set next_dump_any to last+1 so will not influence next

  int writeflag;

  if (ndump && update->restrict_output == 0) {
    for (int idump = 0; idump < ndump; idump++) {
      if (dump[idump]->clearstep || every_dump[idump] == 0)
        modify->clearstep_compute();
      writeflag = 0;
      if (every_dump[idump] && ntimestep % every_dump[idump] == 0 &&
          last_dump[idump] != ntimestep) writeflag = 1;
      if (last_dump[idump] < 0 && dump[idump]->first_flag == 1) writeflag = 1;

      if (writeflag) {
        dump[idump]->write();
        last_dump[idump] = ntimestep;
      }
      if (every_dump[idump])
        next_dump[idump] =
          (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      else {
        bigint nextdump = static_cast<bigint>
          (input->variable->compute_equal(ivar_dump[idump]));
        if (nextdump <= ntimestep)
          error->all(FLERR,"Dump every variable returned a bad timestep");
        next_dump[idump] = nextdump;
      }
      if (dump[idump]->clearstep || every_dump[idump] == 0) {
        if (writeflag) modify->addstep_compute(next_dump[idump]);
        else modify->addstep_compute_all(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  // do not write restart files at start of run
  // set next_restart values to multiple of every or variable value
  // wrap variable eval with clear/add
  // if no restarts, set next_restart to last+1 so will not influence next

  if (restart_flag && update->restrict_output == 0) {
    if (restart_flag_single) {
      if (restart_every_single)
        next_restart_single =
          (ntimestep/restart_every_single)*restart_every_single +
          restart_every_single;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_single = nextrestart;
      }
    } else next_restart_single = update->laststep + 1;
    if (restart_flag_double) {
      if (restart_every_double)
        next_restart_double =
          (ntimestep/restart_every_double)*restart_every_double +
          restart_every_double;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
      }
    } else next_restart_double = update->laststep + 1;
    next_restart = MIN(next_restart_single,next_restart_double);
  } else next_restart = update->laststep + 1;

  // print memory usage unless being called between multiple runs

  if (memflag) memory_usage();

  // set next_thermo to multiple of every or variable eval if var defined
  // insure thermo output on last step of run
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
  // next_dump does not force output on last step of run
  // wrap dumps that invoke computes or eval of variable with clear/add

  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep) {
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->clearstep_compute();
        if (last_dump[idump] != ntimestep) {
          dump[idump]->write();
          last_dump[idump] = ntimestep;
        }
        if (every_dump[idump]) next_dump[idump] += every_dump[idump];
        else {
          bigint nextdump = static_cast<bigint>
            (input->variable->compute_equal(ivar_dump[idump]));
          if (nextdump <= ntimestep)
            error->all(FLERR,"Dump every variable returned a bad timestep");
          next_dump[idump] = nextdump;
        }
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // eval of variable may invoke computes so wrap with clear/add

  if (next_restart == ntimestep) {
    if (next_restart_single == ntimestep) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
      *ptr = '*';
      if (last_restart != ntimestep) restart->write(file);
      delete [] file;
      if (restart_every_single) next_restart_single += restart_every_single;
      else {
        modify->clearstep_compute();
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
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
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
        modify->addstep_compute(next_restart_double);
      }
    }
    last_restart = ntimestep;
    next_restart = MIN(next_restart_single,next_restart_double);
  }

  // insure next_thermo forces output on last step of run
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
   force restart file(s) to be written
   called from PRD and TAD
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_flag_single) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
    *ptr = '*';
    restart->write(file);
    delete [] file;
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
   reset next timestep values for dumps, restart, thermo output
   reset to smallest value >= new timestep
   if next timestep set by variable evaluation,
     eval for ntimestep-1, so current ntimestep can be returned if needed
     no guarantee that variable can be evaluated for ntimestep-1
       if it depends on computes, but live with that rare case for now
------------------------------------------------------------------------- */

void Output::reset_timestep(bigint ntimestep)
{
  next_dump_any = MAXBIGINT;
  for (int idump = 0; idump < ndump; idump++) {
    if (every_dump[idump]) {
      next_dump[idump] = (ntimestep/every_dump[idump])*every_dump[idump];
      if (next_dump[idump] < ntimestep) next_dump[idump] += every_dump[idump];
    } else {
      // ivar_dump may not be initialized
      if (ivar_dump[idump] < 0) {
        ivar_dump[idump] = input->variable->find(var_dump[idump]);
        if (ivar_dump[idump] < 0)
          error->all(FLERR,"Variable name for dump every does not exist");
        if (!input->variable->equalstyle(ivar_dump[idump]))
          error->all(FLERR,"Variable for dump every is invalid style");
      }
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextdump = static_cast<bigint>
        (input->variable->compute_equal(ivar_dump[idump]));
      if (nextdump < ntimestep)
        error->all(FLERR,"Dump every variable returned a bad timestep");
      update->ntimestep++;
      next_dump[idump] = nextdump;
      modify->addstep_compute(next_dump[idump]);
    }
    next_dump_any = MIN(next_dump_any,next_dump[idump]);
  }

  if (restart_flag_single) {
    if (restart_every_single) {
      next_restart_single =
        (ntimestep/restart_every_single)*restart_every_single;
      if (next_restart_single < ntimestep)
        next_restart_single += restart_every_single;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_single));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
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
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_double));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
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
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0)
      error->all(FLERR,"Reuse of dump ID");
  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find dump group ID");
  if (force->inumeric(FLERR,arg[3]) <= 0)
    error->all(FLERR,"Invalid dump frequency");

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    memory->grow(every_dump,max_dump,"output:every_dump");
    memory->grow(next_dump,max_dump,"output:next_dump");
    memory->grow(last_dump,max_dump,"output:last_dump");
    var_dump = (char **)
      memory->srealloc(var_dump,max_dump*sizeof(char *),"output:var_dump");
    memory->grow(ivar_dump,max_dump,"output:ivar_dump");
  }

  // initialize per-dump data to suitable default values

  every_dump[ndump] = 0;
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ivar_dump[ndump] = -1;

  // create the Dump

  if (dump_map->find(arg[2]) != dump_map->end()) {
    DumpCreator dump_creator = (*dump_map)[arg[2]];
    dump[ndump] = dump_creator(lmp, narg, arg);
  }
  else error->all(FLERR,"Unknown dump style");

  every_dump[ndump] = force->inumeric(FLERR,arg[3]);
  if (every_dump[ndump] <= 0) error->all(FLERR,"Illegal dump command");
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ndump++;
}

/* ----------------------------------------------------------------------
   one instance per dump style in style_dump.h
------------------------------------------------------------------------- */

template <typename T>
Dump *Output::dump_creator(LAMMPS *lmp, int narg, char ** arg)
{
  return new T(lmp, narg, arg);
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal dump_modify command");

  // find which dump it is

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Cound not find dump_modify ID");

  dump[idump]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Dump from list of Dumps
------------------------------------------------------------------------- */

void Output::delete_dump(char *id)
{
  // find which dump it is and delete it

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(id,dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Could not find undump ID");

  delete dump[idump];
  delete [] var_dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    every_dump[i-1] = every_dump[i];
    next_dump[i-1] = next_dump[i];
    last_dump[i-1] = last_dump[i];
    var_dump[i-1] = var_dump[i];
    ivar_dump[i-1] = ivar_dump[i];
  }
  ndump--;
}

/* ----------------------------------------------------------------------
   find a dump by ID
   return index of dump or -1 if not found
------------------------------------------------------------------------- */

int Output::find_dump(const char *id)
{
  if (id == NULL) return -1;
  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(id,dump[idump]->id) == 0) break;
  if (idump == ndump) return -1;
  return idump;
}

/* ----------------------------------------------------------------------
   set thermo output frequency from input script
------------------------------------------------------------------------- */

void Output::set_thermo(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal thermo command");

  if (strstr(arg[0],"v_") == arg[0]) {
    delete [] var_thermo;
    int n = strlen(&arg[0][2]) + 1;
    var_thermo = new char[n];
    strcpy(var_thermo,&arg[0][2]);
  } else {
    thermo_every = force->inumeric(FLERR,arg[0]);
    if (thermo_every < 0) error->all(FLERR,"Illegal thermo command");
  }
}

/* ----------------------------------------------------------------------
   new Thermo style
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal thermo_style command");

  // don't allow this so that dipole style can safely allocate inertia vector

  if (domain->box_exist == 0)
    error->all(FLERR,"Thermo_style command before simulation box is defined");

  // warn if previous thermo had been modified via thermo_modify command

  if (thermo->modified && comm->me == 0)
    error->warning(FLERR,"New thermo_style command, "
                   "previous thermo_modify settings will be lost");

  // set thermo = NULL in case new Thermo throws an error

  delete thermo;
  thermo = NULL;
  thermo = new Thermo(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   setup restart capability for single or double output files
   if only one filename and it contains no "*", then append ".*"
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal restart command");

  int every = 0;
  int varflag = 0;

  if (strstr(arg[0],"v_") == arg[0]) varflag = 1;
  else every = force->inumeric(FLERR,arg[0]);

  if (!varflag && every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");

    restart_flag = restart_flag_single = restart_flag_double = 0;
    last_restart = -1;

    delete restart;
    restart = NULL;
    delete [] restart1;
    delete [] restart2a;
    delete [] restart2b;
    restart1 = restart2a = restart2b = NULL;
    delete [] var_restart_single;
    delete [] var_restart_double;
    var_restart_single = var_restart_double = NULL;

    return;
  }

  if (narg < 2) error->all(FLERR,"Illegal restart command");

  int nfile = 0;
  if (narg % 2 == 0) nfile = 1;
  else nfile = 2;

  if (nfile == 1) {
    restart_flag = restart_flag_single = 1;

    if (varflag) {
      delete [] var_restart_single;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_single = new char[n];
      strcpy(var_restart_single,&arg[0][2]);
      restart_every_single = 0;
    } else restart_every_single = every;

    int n = strlen(arg[1]) + 3;
    delete [] restart1;
    restart1 = new char[n];
    strcpy(restart1,arg[1]);
    if (strchr(restart1,'*') == NULL) strcat(restart1,".*");
  }

  if (nfile == 2) {
    restart_flag = restart_flag_double = 1;

    if (varflag) {
      delete [] var_restart_double;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_double = new char[n];
      strcpy(var_restart_double,&arg[0][2]);
      restart_every_double = 0;
    } else restart_every_double = every;

    delete [] restart2a;
    delete [] restart2b;
    restart_toggle = 0;
    int n = strlen(arg[1]) + 3;
    restart2a = new char[n];
    strcpy(restart2a,arg[1]);
    n = strlen(arg[2]) + 1;
    restart2b = new char[n];
    strcpy(restart2b,arg[2]);
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
  if (strstr(arg[1],".mpi")) mpiioflag = 1;
  else mpiioflag = 0;
  if (nfile == 2) {
    if (mpiioflag && !strstr(arg[2],".mpi"))
      error->all(FLERR,"Both restart files must use MPI-IO or neither");
    if (!mpiioflag && strstr(arg[2],".mpi"))
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
  bigint bytes = 0;
  bytes += atom->memory_usage();
  bytes += neighbor->memory_usage();
  bytes += comm->memory_usage();
  bytes += update->memory_usage();
  bytes += force->memory_usage();
  bytes += modify->memory_usage();
  for (int i = 0; i < ndump; i++) bytes += dump[i]->memory_usage();

  double mbytes = bytes/1024.0/1024.0;
  double mbavg,mbmin,mbmax;
  MPI_Reduce(&mbytes,&mbavg,1,MPI_DOUBLE,MPI_SUM,0,world);
  MPI_Reduce(&mbytes,&mbmin,1,MPI_DOUBLE,MPI_MIN,0,world);
  MPI_Reduce(&mbytes,&mbmax,1,MPI_DOUBLE,MPI_MAX,0,world);

  if (comm->me == 0) {
    mbavg /= comm->nprocs;
    if (screen)
      fprintf(screen,"Per MPI rank memory allocation (min/avg/max) = "
              "%.4g | %.4g | %.4g Mbytes\n",mbmin,mbavg,mbmax);
    if (logfile)
      fprintf(logfile,"Per MPI rank memory allocation (min/avg/max) = "
              "%.4g | %.4g | %.4g Mbytes\n",mbmin,mbavg,mbmax);
  }
}
