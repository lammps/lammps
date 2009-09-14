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
#include "atom.h"
#include "neighbor.h"
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

#define DumpInclude
#include "style.h"
#undef DumpInclude

using namespace LAMMPS_NS;

#define DELTA 1

#define MYMIN(a,b) ((a) < (b) ? (a) : (b))
#define MYMAX(a,b) ((a) > (b) ? (a) : (b))

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

  ndump = 0;
  max_dump = 0;
  next_dump = NULL;
  last_dump = NULL;
  dump_every = NULL;
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

  memory->sfree(next_dump);
  memory->sfree(last_dump);
  memory->sfree(dump_every);
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
  for (int i = 0; i < ndump; i++) dump[i]->init();
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do thermo last, so will print after memory_usage
------------------------------------------------------------------------- */

void Output::setup(int flag)
{
  int ntimestep = update->ntimestep;

  // perform dump at start of run if current timestep is multiple of every
  //   and last dump was not on this timestep
  // set next_dump to multiple of every
  // will not write on last step of run unless multiple of every
  // set next_dump_any to smallest next_dump
  // if no dumps, set next_dump_any to last+1 so will not influence next
  // dump custom/cfg may invoke computes so wrap with clear/add

  if (ndump && update->restrict_output == 0) {
    for (int idump = 0; idump < ndump; idump++) {
      if (strcmp(dump[idump]->style,"custom") == 0 ||
	  strcmp(dump[idump]->style,"cfg") == 0)
	modify->clearstep_compute();
      if ((ntimestep % dump_every[idump] == 0 && 
	  last_dump[idump] != ntimestep) || last_dump[idump] < 0) {
	dump[idump]->write();
	last_dump[idump] = ntimestep;
      }
      next_dump[idump] = 
	(ntimestep/dump_every[idump])*dump_every[idump] + dump_every[idump];
      if (strcmp(dump[idump]->style,"custom") == 0 ||
	  strcmp(dump[idump]->style,"cfg") == 0)
	modify->addstep_compute(next_dump[idump]);
      if (idump) next_dump_any = MYMIN(next_dump_any,next_dump[idump]);
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
    next_thermo = MYMIN(next_thermo,update->laststep);
  } else next_thermo = update->laststep;

  modify->addstep_compute(next_thermo);

  // next = next timestep any output will be done

  next = MYMIN(next_dump_any,next_restart);
  next = MYMIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last doesn't
   do dump/restart before thermo so thermo CPU time will include them
------------------------------------------------------------------------- */

void Output::write(int ntimestep)
{
  // next_dump does not force output on last step of run
  // dump custom/cfg may invoke computes so wrap with clear/add

  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep && last_dump[idump] != ntimestep) {
        if (strcmp(dump[idump]->style,"custom") == 0 ||
	    strcmp(dump[idump]->style,"cfg") == 0)
          modify->clearstep_compute();
	dump[idump]->write();
	last_dump[idump] = ntimestep;
	next_dump[idump] += dump_every[idump];
        if (strcmp(dump[idump]->style,"custom") == 0 ||
	    strcmp(dump[idump]->style,"cfg") == 0)
	  modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MYMIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename

  if (next_restart == ntimestep && last_restart != ntimestep) {
    if (restart_toggle == 0) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s%d%s",restart1,ntimestep,ptr+1);
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
    next_thermo += thermo_every;
    next_thermo = MYMIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  }

  // next = next timestep any output will be done

  next = MYMIN(next_dump_any,next_restart);
  next = MYMIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   force a snapshot to be written for all dumps
------------------------------------------------------------------------- */

void Output::write_dump(int ntimestep)
{
  for (int idump = 0; idump < ndump; idump++) {
    dump[idump]->write();
    last_dump[idump] = ntimestep;
  }
}

/* ----------------------------------------------------------------------
   force a restart file to be written
------------------------------------------------------------------------- */

void Output::write_restart(int ntimestep)
{
  if (restart_toggle == 0) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s%d%s",restart1,ntimestep,ptr+1);
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
  if (narg < 5) error->all("Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) error->all("Reuse of dump ID");
  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find dump group ID");
  if (atoi(arg[3]) <= 0) error->all("Invalid dump frequency");

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    dump_every = (int *)
      memory->srealloc(dump_every,max_dump*sizeof(int *),"output:dump_every");
    next_dump = (int *)
      memory->srealloc(next_dump,max_dump*sizeof(int *),"output:next_dump");
    last_dump = (int *)
      memory->srealloc(last_dump,max_dump*sizeof(int *),"output:last_dump");
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DumpClass
#define DumpStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) dump[ndump] = new Class(lmp,narg,arg);
#include "style.h"
#undef DumpClass

  else error->all("Invalid dump style");

  dump_every[ndump] = atoi(arg[3]);
  if (dump_every[ndump] <= 0) error->all("Illegal dump command");
  last_dump[ndump] = -1;
  ndump++;
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump 
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal dump_modify command");

  // find which dump it is

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) break;
  if (idump == ndump) error->all("Cound not find dump_modify ID");

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
  if (idump == ndump) error->all("Could not find undump ID");

  delete dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    dump_every[i-1] = dump_every[i];
    next_dump[i-1] = next_dump[i];
    last_dump[i-1] = last_dump[i];
  }
  ndump--;
}

/* ----------------------------------------------------------------------
   new Thermo style 
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal thermo_style command");

  // don't allow this so that dipole style can safely allocate inertia vector

  if (domain->box_exist == 0) 
    error->all("Thermo_style command before simulation box is defined");

  // warn if previous thermo had been modified via thermo_modify command

  if (thermo->modified && comm->me == 0)
    error->warning("New thermo_style command, "
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
  if (narg < 1) error->all("Illegal restart command");

  if (restart) delete restart;
  delete [] restart1;
  delete [] restart2;
  restart = NULL;
  restart1 = restart2 = NULL;
  last_restart = -1;

  restart_every = atoi(arg[0]);
  if (restart_every == 0) {
    if (narg != 1) error->all("Illegal restart command");
    return;
  }

  restart = new WriteRestart(lmp);

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
  } else error->all("Illegal restart command");
}

/* ----------------------------------------------------------------------
   sum and print memory usage
   is only memory on proc 0, not averaged across procs
------------------------------------------------------------------------- */

void Output::memory_usage()
{
  double bytes = 0.0;
  bytes += atom->memory_usage();
  bytes += neighbor->memory_usage();
  bytes += comm->memory_usage();
  bytes += update->memory_usage();
  bytes += force->memory_usage();
  bytes += modify->memory_usage();
  for (int i = 0; i < ndump; i++) bytes += dump[i]->memory_usage();

  double mbytes = bytes/1024.0/1024.0;

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Memory usage per processor = %g Mbytes\n",mbytes);
    if (logfile) 
      fprintf(logfile,"Memory usage per processor = %g Mbytes\n",mbytes);
  }
}
