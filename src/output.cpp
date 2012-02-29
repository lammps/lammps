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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
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
#include "accelerator_cuda.h"
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
  modify->add_compute(3,newarg,lmp->suffix);

  newarg[0] = (char *) "thermo_press";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "thermo_temp";
  modify->add_compute(4,newarg,lmp->suffix);

  newarg[0] = (char *) "thermo_pe";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg,lmp->suffix);

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

  restart = NULL;
  restart1 = restart2 = NULL;
  restart_every = 0;
  last_restart = -1;
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

  delete restart;
  delete [] restart1;
  delete [] restart2;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  thermo->init();
  if (thermo_every) delete [] var_thermo;
  else if (var_thermo) {
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
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do thermo last, so will print after memory_usage
------------------------------------------------------------------------- */

void Output::setup(int flag)
{
  bigint ntimestep = update->ntimestep;

  // perform dump at start of run if current timestep is multiple of every
  //   and last dump was not on this timestep
  // set next_dump to multiple of every
  // will not write on last step of run unless multiple of every
  // set next_dump_any to smallest next_dump
  // if no dumps, set next_dump_any to last+1 so will not influence next
  // wrap dumps that invoke computes with clear/add
  // if dump not written now, add_all on future step since clear/add is noop

  int writeflag;

  if (ndump && update->restrict_output == 0) {
    for (int idump = 0; idump < ndump; idump++) {
      if (dump[idump]->clearstep) modify->clearstep_compute();
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
	int nextdump = static_cast<int> 
	  (input->variable->compute_equal(ivar_dump[idump]));
	if (nextdump <= ntimestep)
	  error->all(FLERR,"Dump every variable returned a bad timestep");
	next_dump[idump] = nextdump;
      }
      if (dump[idump]->clearstep) {
	if (writeflag) modify->addstep_compute(next_dump[idump]);
	else modify->addstep_compute_all(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  // do not write a restart file at start of run
  // set next_restart to multiple of every
  // will not write on last step of run unless multiple of every
  // if every = 0, set next_restart to last+1 so will not influence next

  if (restart_every && update->restrict_output == 0)
    next_restart = (ntimestep/restart_every)*restart_every + restart_every;
  else next_restart = update->laststep + 1;

  // print memory usage unless being called between multiple runs

  if (flag) memory_usage();

  // always do thermo with header at start of run
  // set next_thermo to multiple of every or last step of run (if smaller)
  // if every = 0, set next_thermo to last step of run
  // thermo may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  thermo->header();
  thermo->compute(0);
  last_thermo = ntimestep;

  if (thermo_every) {
    next_thermo = (ntimestep/thermo_every)*thermo_every + thermo_every;
    next_thermo = MIN(next_thermo,update->laststep);
  } else if (var_thermo) {
    next_thermo = static_cast<int> 
      (input->variable->compute_equal(ivar_thermo));
    if (next_thermo <= ntimestep)
      error->all(FLERR,"Thermo every variable returned a bad timestep");
  } else next_thermo = update->laststep;

  modify->addstep_compute(next_thermo);

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last doesn't
   do dump/restart before thermo so thermo CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  // next_dump does not force output on last step of run
  // wrap dumps that invoke computes with clear/add
  // download data from GPU if necessary

  if (next_dump_any == ntimestep) {

    if (lmp->cuda && !lmp->cuda->oncpu) lmp->cuda->downloadAll();    
    
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep && last_dump[idump] != ntimestep) {
        if (dump[idump]->clearstep) modify->clearstep_compute();
	dump[idump]->write();
	last_dump[idump] = ntimestep;
	if (every_dump[idump]) next_dump[idump] += every_dump[idump];
	else {
	  int nextdump = static_cast<int> 
	    (input->variable->compute_equal(ivar_dump[idump]));
	  if (nextdump <= ntimestep)
	    error->all(FLERR,"Dump every variable returned a bad timestep");
	  next_dump[idump] = nextdump;
	}
        if (dump[idump]->clearstep) modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // download data from GPU if necessary

  if (next_restart == ntimestep && last_restart != ntimestep) {

    if (lmp->cuda && !lmp->cuda->oncpu) lmp->cuda->downloadAll();    
    
    if (restart_toggle == 0) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
      *ptr = '*';
      restart->write(file);
      delete [] file;
    } else if (restart_toggle == 1) {
      restart->write(restart1);
      restart_toggle = 2;
    } else if (restart_toggle == 2) {
      restart->write(restart2);
      restart_toggle = 1;
    }
    last_restart = ntimestep;
    next_restart += restart_every;
  }

  // insure next_thermo forces output on last step of run
  // thermo may invoke computes so wrap with clear/add

  if (next_thermo == ntimestep && last_thermo != ntimestep) {
    modify->clearstep_compute();
    thermo->compute(1);
    last_thermo = ntimestep;
    if (thermo_every) next_thermo += thermo_every;
    else if (var_thermo) {
      next_thermo = static_cast<int> 
	(input->variable->compute_equal(ivar_thermo));
      if (next_thermo <= ntimestep)
	error->all(FLERR,"Thermo every variable returned a bad timestep");
    } else next_thermo = update->laststep;
    next_thermo = MIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  }

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   force a snapshot to be written for all dumps
------------------------------------------------------------------------- */

void Output::write_dump(bigint ntimestep)
{
  for (int idump = 0; idump < ndump; idump++) {
    dump[idump]->write();
    last_dump[idump] = ntimestep;
  }
}

/* ----------------------------------------------------------------------
   force a restart file to be written
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_toggle == 0) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
    *ptr = '*';
    restart->write(file);
    delete [] file;
  } else if (restart_toggle == 1) {
    restart->write(restart1);
    restart_toggle = 2;
  } else if (restart_toggle == 2) {
    restart->write(restart2);
    restart_toggle = 1;
  }

  last_restart = ntimestep;
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
  if (atoi(arg[3]) <= 0) error->all(FLERR,"Invalid dump frequency");

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

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) dump[ndump] = new Class(lmp,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  every_dump[ndump] = atoi(arg[3]);
  if (every_dump[ndump] <= 0) error->all(FLERR,"Illegal dump command");
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ndump++;
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
   setup restart capability
   if only one filename and it contains no "*", then append ".*"
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal restart command");

  if (restart) delete restart;
  delete [] restart1;
  delete [] restart2;
  restart = NULL;
  restart1 = restart2 = NULL;
  last_restart = -1;

  restart_every = atoi(arg[0]);
  if (restart_every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");
    return;
  }

  restart = new WriteRestart(lmp);

  if (narg != 2 && narg != 3) error->all(FLERR,"Illegal restart command");

  int n = strlen(arg[1]) + 3;
  restart1 = new char[n];
  strcpy(restart1,arg[1]);

  if (narg == 2) {
    restart_toggle = 0;
    restart2 = NULL;
    if (strchr(restart1,'*') == NULL) strcat(restart1,".*");
  } else if (narg == 3) {
    restart_toggle = 1;
    n = strlen(arg[2]) + 1;
    restart2 = new char[n];
    strcpy(restart2,arg[2]);
  } else error->all(FLERR,"Illegal restart command");
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
  for (int i = 0; i < ndump; i++) dump[i]->memory_usage();

  double mbytes = bytes/1024.0/1024.0;

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Memory usage per processor = %g Mbytes\n",mbytes);
    if (logfile) 
      fprintf(logfile,"Memory usage per processor = %g Mbytes\n",mbytes);
  }
}
