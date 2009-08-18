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
#include "string.h"
#include "thermo.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// customize a new keyword by adding to this list:

// step, atoms, cpu, temp, press, pe, ke, etotal, enthalpy
// evdwl, ecoul, epair, ebond, eangle, edihed, eimp, emol, elong, etail
// vol, lx, ly, lz, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
// pxx, pyy, pzz, pxy, pxz, pyz

// customize a new thermo style by adding a DEFINE to this list

#define ONE "step temp epair emol etotal press"
#define MULTI "etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press"

enum{IGNORE,WARN,ERROR};           // same as write_restart.cpp
enum{ONELINE,MULTILINE};
enum{INT,FLOAT};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

#define MAXLINE 1024
#define DELTA 8

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Thermo::Thermo(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  // set thermo_modify defaults

  modified = 0;
  normuserflag = 0;
  lineflag = ONELINE;
  lostflag = ERROR;
  lostbefore = 0;
  flushflag = 0;

  // set style and corresponding lineflag
  // custom style builds its own line of keywords
  // customize a new thermo style by adding to if statement

  line = new char[MAXLINE];

  if (strcmp(style,"one") == 0) {
    strcpy(line,ONE);
  } else if (strcmp(style,"multi") == 0) {
    strcpy(line,MULTI);
    lineflag = MULTILINE;

  } else if (strcmp(style,"custom") == 0) {
    if (narg == 1) error->all("Illegal thermo style custom command");
    line[0] = '\0';
    for (int iarg = 1; iarg < narg; iarg++) {
      strcat(line,arg[iarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';

  } else error->all("Illegal thermo style command");

  // ptrs, flags, IDs for compute objects thermo may use or create

  temperature = NULL;
  pressure = NULL;
  pe = NULL;

  index_temp = index_press_scalar = index_press_vector = index_pe = -1;

  id_temp = (char *) "thermo_temp";
  id_press = (char *) "thermo_press";
  id_pe = (char *) "thermo_pe";

  // count fields in line
  // allocate per-field memory
  // process line of keywords

  nfield_initial = atom->count_words(line);
  allocate();
  parse_fields(line);

  // format strings

  format_multi = (char *) "---------------- Step %8d ----- "
                          "CPU = %11.4f (sec) ----------------";
  format_int_one_def = (char *) "%8d";
  format_int_multi_def = (char *) "%14d";
  format_g_def = (char *) "%12.8g";
  format_f_def = (char *) "%14.4f";
  format_int_user = NULL;
  format_float_user = NULL;
}

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
  delete [] style;
  delete [] line;

  deallocate();

  // format strings

  delete [] format_int_user;
  delete [] format_float_user;
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
  int i,n;

  // set normvalue to default setting unless user has specified it

  if (normuserflag) normvalue = normuser;
  else if (strcmp(update->unit_style,"lj") == 0) normvalue = 1;
  else normvalue = 0;

  // add Volume field if volume changes and not style = custom
  // this check must come after domain init, so box_change is set

  nfield = nfield_initial;
  if (domain->box_change && strcmp(style,"custom") != 0)
    addfield("Volume",&Thermo::compute_vol,FLOAT);

  // set format string for each field
  // include keyword if lineflag = MULTILINE
  // add '/n' every 3 values if lineflag = MULTILINE
  // add trailing '/n' to last value

  char *ptr;
  for (i = 0; i < nfield; i++) {
    format[i][0] = '\0';
    if (lineflag == MULTILINE && i % 3 == 0) strcat(format[i],"\n");

    if (format_user[i]) ptr = format_user[i];
    else if (vtype[i] == INT && format_int_user) ptr = format_int_user;
    else if (vtype[i] == INT && lineflag == ONELINE) ptr = format_int_one_def;
    else if (vtype[i] == INT && lineflag == MULTILINE) 
      ptr = format_int_multi_def;
    else if (vtype[i] == FLOAT && format_float_user) ptr = format_float_user;
    else if (lineflag == ONELINE) ptr = format_g_def;
    else if (lineflag == MULTILINE) ptr = format_f_def;

    n = strlen(format[i]);
    if (lineflag == ONELINE) sprintf(&format[i][n],"%s ",ptr);
    else sprintf(&format[i][n],"%-8s = %s ",keyword[i],ptr);

    if (i == nfield-1) strcat(format[i],"\n");
  }

  // find current ptr for each Compute ID

  int icompute;
  for (i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all("Could not find thermo compute ID");
    computes[i] = modify->compute[icompute];
  }

  // find current ptr for each Fix ID
  // check that fix frequency is acceptable with thermo output frequency

  int ifix;
  for (i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all("Could not find thermo fix ID");
    fixes[i] = modify->fix[ifix];
    if (output->thermo_every % fixes[i]->scalar_vector_freq)
      error->all("Thermo and fix not computed at compatible times");
  }

  // find current ptr for each Variable ID

  int ivariable;
  for (i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) 
      error->all("Could not find thermo custom variable name");
    variables[i] = ivariable;
  }

  // set ptrs to keyword-specific Compute objects

  if (index_temp >= 0) temperature = computes[index_temp];
  if (index_press_scalar >= 0) pressure = computes[index_press_scalar];
  if (index_press_vector >= 0) pressure = computes[index_press_vector];
  if (index_pe >= 0) pe = computes[index_pe];
}

/* ---------------------------------------------------------------------- */

void Thermo::header()
{
  if (lineflag == MULTILINE) return;

  int loc = 0;
  for (int i = 0; i < nfield; i++)
    loc += sprintf(&line[loc],"%s ",keyword[i]);
  sprintf(&line[loc],"\n");
  
  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) fprintf(logfile,line);
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute(int flag)
{
  int i;

  firststep = flag;
  int ntimestep = update->ntimestep;

  // check for lost atoms
  // turn off normflag if natoms = 0 to avoid divide by 0

  natoms = lost_check();
  if (natoms == 0) normflag = 0;
  else normflag = normvalue;

  // invoke Compute methods needed for thermo keywords
  // which = 0 is global scalar, which = 1 is global vector

  for (i = 0; i < ncompute; i++)
    if (compute_which[i] == 0) {
      if (!(computes[i]->invoked_flag & INVOKED_SCALAR)) {
	computes[i]->compute_scalar();
	computes[i]->invoked_flag |= INVOKED_SCALAR;
      }
    } else {
      if (!(computes[i]->invoked_flag & INVOKED_VECTOR)) {
	computes[i]->compute_vector();
	computes[i]->invoked_flag |= INVOKED_VECTOR;
      }
    }

  // if lineflag = MULTILINE, prepend step/cpu header line

  int loc = 0;
  if (lineflag == MULTILINE) {
    double cpu;
    if (flag) cpu = timer->elapsed(TIME_LOOP);
    else cpu = 0.0;
    loc = sprintf(&line[loc],format_multi,ntimestep,cpu);
  }

  // add each thermo value to line with its specific format

  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == INT) loc += sprintf(&line[loc],format[ifield],ivalue);
    else loc += sprintf(&line[loc],format[ifield],dvalue);
  }

  // kludge for RedStorm timing issue
  // if (ntimestep == 100) return;

  // print line to screen and logfile

  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) {
      fprintf(logfile,line);
      if (flushflag) fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   check for lost atoms, return current number of atoms
------------------------------------------------------------------------- */

double Thermo::lost_check()
{
  // ntotal = current # of atoms

  double ntotal;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);
  if (ntotal == atom->natoms) return ntotal;

  // if not checking or already warned, just return

  if (lostflag == IGNORE) return ntotal;
  if (lostflag == WARN && lostbefore == 1) return ntotal;

  // error message

  if (lostflag == ERROR) {
    char str[128];
    sprintf(str,"Lost atoms: original %.15g current %.15g",
	    atom->natoms,ntotal);
    error->all(str);
  }

  // warning message

  char str[128];
  sprintf(str,"Lost atoms: original %.15g current %.15g",atom->natoms,ntotal);
  if (me == 0) error->warning(str);
  lostbefore = 1;
  return ntotal;
}

/* ----------------------------------------------------------------------
   modify thermo parameters
------------------------------------------------------------------------- */

void Thermo::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal thermo_modify command");

  modified = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_temp < 0) error->all("Thermo style does not use temp");
      delete [] id_compute[index_temp];
      int n = strlen(arg[iarg+1]) + 1;
      id_compute[index_temp] = new char[n];
      strcpy(id_compute[index_temp],arg[iarg+1]);

      int icompute = modify->find_compute(arg[iarg+1]);
      if (icompute < 0) 
	error->all("Could not find thermo_modify temperature ID");
      temperature = modify->compute[icompute];

      if (temperature->tempflag == 0)
	error->all("Thermo_modify temperature ID does not "
		   "compute temperature");
      if (temperature->igroup != 0 && comm->me == 0)
	error->warning("Temperature for thermo pressure is not for group all");

      // reset id_temp of pressure to new temperature ID
      // either pressure currently being used by thermo or "thermo_press"

      if (index_press_scalar >= 0) {
	icompute = modify->find_compute(id_compute[index_press_scalar]);
	if (icompute < 0) error->all("Pressure ID for thermo does not exist");
      } else if (index_press_vector >= 0) {
	icompute = modify->find_compute(id_compute[index_press_vector]);
	if (icompute < 0) error->all("Pressure ID for thermo does not exist");
      } else icompute = modify->find_compute((char *) "thermo_press");

      modify->compute[icompute]->reset_extra_compute_fix(arg[iarg+1]);

      iarg += 2;

    } else if (strcmp(arg[iarg],"press") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (index_press_scalar < 0 && index_press_vector < 0)
	error->all("Thermo style does not use press");

      if (index_press_scalar >= 0) {
	delete [] id_compute[index_press_scalar];
	int n = strlen(arg[iarg+1]) + 1;
	id_compute[index_press_scalar] = new char[n];
	strcpy(id_compute[index_press_scalar],arg[iarg+1]);
      }
      if (index_press_vector >= 0) {
	delete [] id_compute[index_press_vector];
	int n = strlen(arg[iarg+1]) + 1;
	id_compute[index_press_vector] = new char[n];
	strcpy(id_compute[index_press_vector],arg[iarg+1]);
      }

      int icompute = modify->find_compute(arg[iarg+1]);
      if (icompute < 0) error->all("Could not find thermo_modify pressure ID");
      pressure = modify->compute[icompute];

      if (pressure->pressflag == 0)
	error->all("Thermo_modify pressure ID does not compute pressure");

      iarg += 2;

    } else if (strcmp(arg[iarg],"lost") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"ignore") == 0) lostflag = IGNORE;
      else if (strcmp(arg[iarg+1],"warn") == 0) lostflag = WARN;
      else if (strcmp(arg[iarg+1],"error") == 0) lostflag = ERROR;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      normuserflag = 1;
      if (strcmp(arg[iarg+1],"no") == 0) normuser = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) normuser = 1;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) flushflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) flushflag = 1;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"line") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"one") == 0) lineflag = ONELINE;
      else if (strcmp(arg[iarg+1],"multi") == 0) lineflag = MULTILINE;
      else error->all("Illegal thermo_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+3 > narg) error->all("Illegal thermo_modify command");
      if (strcmp(arg[iarg+1],"int") == 0) {
	if (format_int_user) delete [] format_int_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_int_user = new char[n];
	strcpy(format_int_user,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"float") == 0) {
	if (format_float_user) delete [] format_float_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_float_user = new char[n];
	strcpy(format_float_user,arg[iarg+2]);
      } else {
	int i = atoi(arg[iarg+1]) - 1;
	if (i < 0 || i >= nfield_initial)
	  error->all("Illegal thermo_modify command");
	if (format_user[i]) delete [] format_user[i];
	int n = strlen(arg[iarg+2]) + 1;
	format_user[i] = new char[n];
	strcpy(format_user[i],arg[iarg+2]);
      }
      iarg += 3;

    } else error->all("Illegal thermo_modify command");
  }
}

/* ----------------------------------------------------------------------
   allocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::allocate()
{
  // n = specified fields + Volume field (added at run time)

  int n = nfield_initial + 1;

  keyword = new char*[n];
  for (int i = 0; i < n; i++) keyword[i] = new char[32];
  vfunc = new FnPtr[n];
  vtype = new int[n];

  format = new char*[n];
  for (int i = 0; i < n; i++) format[i] = new char[32];
  format_user = new char*[n];
  for (int i = 0; i < n; i++) format_user[i] = NULL;

  field2index = new int[n];
  argindex = new int[n];

  // factor of 3 is max number of computes a single field can add

  ncompute = 0;
  id_compute = new char*[3*n];
  compute_which = new int[3*n];
  computes = new Compute*[3*n];

  nfix = 0;
  id_fix = new char*[n];
  fixes = new Fix*[n];

  nvariable = 0;
  id_variable = new char*[n];
  variables = new int[n];
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::deallocate()
{
  int n = nfield_initial + 1;

  for (int i = 0; i < n; i++) delete [] keyword[i];
  delete [] keyword;
  delete [] vfunc;
  delete [] vtype;

  for (int i = 0; i < n; i++) delete [] format[i];
  delete [] format;
  for (int i = 0; i < n; i++) delete [] format_user[i];
  delete [] format_user;

  delete [] field2index;
  delete [] argindex;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  delete [] id_compute;
  delete [] compute_which;
  delete [] computes;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  delete [] id_fix;
  delete [] fixes;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  delete [] id_variable;
  delete [] variables;
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, pe, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(char *str)
{
  nfield = 0;

  // customize a new keyword by adding to if statement

  char *word = strtok(str," \0");
  while (word) {

    if (strcmp(word,"step") == 0) {
      addfield("Step",&Thermo::compute_step,INT);
    } else if (strcmp(word,"atoms") == 0) {
      addfield("Atoms",&Thermo::compute_atoms,INT);
    } else if (strcmp(word,"cpu") == 0) {
      addfield("CPU",&Thermo::compute_cpu,FLOAT);

    } else if (strcmp(word,"temp") == 0) {
      addfield("Temp",&Thermo::compute_temp,FLOAT);
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"press") == 0) {
      addfield("Press",&Thermo::compute_press,FLOAT);
      index_press_scalar = add_compute(id_press,0);
    } else if (strcmp(word,"pe") == 0) {
      addfield("PotEng",&Thermo::compute_pe,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"ke") == 0) {
      addfield("KinEng",&Thermo::compute_ke,FLOAT);
      index_temp = add_compute(id_temp,0);
    } else if (strcmp(word,"etotal") == 0) {
      addfield("TotEng",&Thermo::compute_etotal,FLOAT);
      index_temp = add_compute(id_temp,0);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"enthalpy") == 0) {
      addfield("Enthalpy",&Thermo::compute_enthalpy,FLOAT);
      index_temp = add_compute(id_temp,0);
      index_press_scalar = add_compute(id_press,0);
      index_pe = add_compute(id_pe,0);

    } else if (strcmp(word,"evdwl") == 0) {
      addfield("E_vdwl",&Thermo::compute_evdwl,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"ecoul") == 0) {
      addfield("E_coul",&Thermo::compute_ecoul,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"epair") == 0) {
      addfield("E_pair",&Thermo::compute_epair,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"ebond") == 0) {
      addfield("E_bond",&Thermo::compute_ebond,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"eangle") == 0) {
      addfield("E_angle",&Thermo::compute_eangle,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"edihed") == 0) {
      addfield("E_dihed",&Thermo::compute_edihed,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"eimp") == 0) {
      addfield("E_impro",&Thermo::compute_eimp,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"emol") == 0) {
      addfield("E_mol",&Thermo::compute_emol,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"elong") == 0) {
      addfield("E_long",&Thermo::compute_elong,FLOAT);
      index_pe = add_compute(id_pe,0);
    } else if (strcmp(word,"etail") == 0) {
      addfield("E_tail",&Thermo::compute_etail,FLOAT);
      index_pe = add_compute(id_pe,0);

    } else if (strcmp(word,"vol") == 0) {
      addfield("Volume",&Thermo::compute_vol,FLOAT);
    } else if (strcmp(word,"lx") == 0) {
      addfield("Lx",&Thermo::compute_lx,FLOAT);
    } else if (strcmp(word,"ly") == 0) {
      addfield("Ly",&Thermo::compute_ly,FLOAT);
    } else if (strcmp(word,"lz") == 0) {
      addfield("Lz",&Thermo::compute_lz,FLOAT);

    } else if (strcmp(word,"xlo") == 0) {
      addfield("Xlo",&Thermo::compute_xlo,FLOAT);
    } else if (strcmp(word,"xhi") == 0) {
      addfield("Xhi",&Thermo::compute_xhi,FLOAT);
    } else if (strcmp(word,"ylo") == 0) {
      addfield("Ylo",&Thermo::compute_ylo,FLOAT);
    } else if (strcmp(word,"yhi") == 0) {
      addfield("Yhi",&Thermo::compute_yhi,FLOAT);
    } else if (strcmp(word,"zlo") == 0) {
      addfield("Zlo",&Thermo::compute_zlo,FLOAT);
    } else if (strcmp(word,"zhi") == 0) {
      addfield("Zhi",&Thermo::compute_zhi,FLOAT);

    } else if (strcmp(word,"xy") == 0) {
      addfield("Xy",&Thermo::compute_xy,FLOAT);
    } else if (strcmp(word,"xz") == 0) {
      addfield("Xz",&Thermo::compute_xz,FLOAT);
    } else if (strcmp(word,"yz") == 0) {
      addfield("Yz",&Thermo::compute_yz,FLOAT);

    } else if (strcmp(word,"pxx") == 0) {
      addfield("Pxx",&Thermo::compute_pxx,FLOAT);
      index_press_vector = add_compute(id_press,1);
    } else if (strcmp(word,"pyy") == 0) {
      addfield("Pyy",&Thermo::compute_pyy,FLOAT);
      index_press_vector = add_compute(id_press,1);
    } else if (strcmp(word,"pzz") == 0) {
      addfield("Pzz",&Thermo::compute_pzz,FLOAT);
      index_press_vector = add_compute(id_press,1);
    } else if (strcmp(word,"pxy") == 0) {
      addfield("Pxy",&Thermo::compute_pxy,FLOAT);
      index_press_vector = add_compute(id_press,1);
    } else if (strcmp(word,"pxz") == 0) {
      addfield("Pxz",&Thermo::compute_pxz,FLOAT);
      index_press_vector = add_compute(id_press,1);
    } else if (strcmp(word,"pyz") == 0) {
      addfield("Pyz",&Thermo::compute_pyz,FLOAT);
      index_press_vector = add_compute(id_press,1);

    // compute value = c_ID, fix value = f_ID, variable value = v_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // copy = at most 8 chars of ID to pass to addfield

    } else if ((strncmp(word,"c_",2) == 0) || (strncmp(word,"f_",2) == 0) ||
	       (strncmp(word,"v_",2) == 0)) {

      int n = strlen(word);
      char *id = new char[n];
      strcpy(id,&word[2]);
      char copy[9];
      strncpy(copy,id,8);
      copy[8] = '\0';

      char *ptr = strchr(id,'[');
      if (ptr) {
	if (id[strlen(id)-1] != ']')
	  error->all("Invalid keyword in thermo_style custom command");
	argindex[nfield] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nfield] = 0;

      if (word[0] == 'c') {
	n = modify->find_compute(id);
	if (n < 0) error->all("Could not find thermo custom compute ID");
	if (argindex[nfield] == 0 && modify->compute[n]->scalar_flag == 0)
	  error->all("Thermo compute ID does not compute scalar info");
	if (argindex[nfield] > 0 && modify->compute[n]->vector_flag == 0)
	  error->all("Thermo compute ID does not compute vector info");
	if (argindex[nfield] > 0 && 
	    argindex[nfield] > modify->compute[n]->size_vector)
	  error->all("Thermo compute ID vector is not large enough");

	field2index[nfield] = add_compute(id,MIN(argindex[nfield],1));
	addfield(copy,&Thermo::compute_compute,FLOAT);

      } else if (word[0] == 'f') {
	n = modify->find_fix(id);
	if (n < 0) error->all("Could not find thermo custom fix ID");
	if (argindex[nfield] == 0 && modify->fix[n]->scalar_flag == 0)
	  error->all("Thermo fix ID does not compute scalar info");
	if (argindex[nfield] > 0 && modify->fix[n]->vector_flag == 0)
	  error->all("Thermo fix ID does not compute vector info");
	if (argindex[nfield] > 0 && 
	    argindex[nfield] > modify->fix[n]->size_vector)
	  error->all("Thermo fix ID vector is not large enough");

	field2index[nfield] = add_fix(id);
	addfield(copy,&Thermo::compute_fix,FLOAT);

      } else if (word[0] == 'v') {
	n = input->variable->find(id);
	if (n < 0) error->all("Could not find thermo custom variable name");
	if (input->variable->equalstyle(n) == 0)
	  error->all("Thermo custom variable is not equal-style variable");

	field2index[nfield] = add_variable(id);
	addfield(copy,&Thermo::compute_variable,FLOAT);
      }

      delete [] id;

    } else error->all("Invalid keyword in thermo_style custom command");

    word = strtok(NULL," \0");
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Thermo::addfield(const char *key, FnPtr func, int typeflag)
{
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   add compute ID to list of Compute objects to call
   return location of where this Compute is in list
   if already in list with same which, do not add, just return index
------------------------------------------------------------------------- */

int Thermo::add_compute(const char *id, int which)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if ((strcmp(id,id_compute[icompute]) == 0) && 
	which == compute_which[icompute]) break;
  if (icompute < ncompute) return icompute;

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  compute_which[ncompute] = which;
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add fix ID to list of Fix objects to call
------------------------------------------------------------------------- */

int Thermo::add_fix(const char *id)
{
  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add variable ID to list of Variables to evaluate
------------------------------------------------------------------------- */

int Thermo::add_variable(const char *id)
{
  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   compute a single thermodyanmic value, word is any keyword in custom list
   called when a variable is evaluated by Variable class
   return value as double in answer
   return 0 if str is recoginzed keyword, 1 if unrecognized
   customize a new keyword by adding to if statement
------------------------------------------------------------------------- */

int Thermo::evaluate_keyword(char *word, double *answer)
{
  // invoke a lo-level thermo routine to compute the variable value
  // if keyword requires a compute, error if thermo doesn't use the compute
  // if inbetween runs and needed compute is not current, error
  // if in middle of run and needed compute is not current, invoke it
  // for keywords that use pe indirectly (evdwl, ebond, etc):
  //   check if energy was tallied on this timestep and set pe->invoked_flag 
  //   this will trigger next timestep for energy tallying via addstep()

  if (strcmp(word,"step") == 0) {
    compute_step();
    dvalue = ivalue;

  } else if (strcmp(word,"atoms") == 0) {
    compute_atoms();
    dvalue = ivalue;

  } else if (strcmp(word,"cpu") == 0) {
    if (update->whichflag == 0) firststep = 0;
    compute_cpu();

  } else if (strcmp(word,"temp") == 0) {
    if (!temperature)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init temp");	  
    if (update->whichflag == 0) {
      if (temperature->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(temperature->invoked_flag & INVOKED_SCALAR)) {
      temperature->compute_scalar();
      temperature->invoked_flag |= INVOKED_SCALAR;
    }
    compute_temp();
    
  } else if (strcmp(word,"press") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_SCALAR)) {
      pressure->compute_scalar();
      pressure->invoked_flag |= INVOKED_SCALAR;
    }
    compute_press();

  } else if (strcmp(word,"pe") == 0) {
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    if (update->whichflag == 0) {
      if (pe->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pe->invoked_flag & INVOKED_SCALAR)) {
      pe->compute_scalar();
      pe->invoked_flag |= INVOKED_SCALAR;
    }
    compute_pe();

  } else if (strcmp(word,"ke") == 0) {
    if (!temperature)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init temp");
    if (update->whichflag == 0) {
      if (temperature->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(temperature->invoked_flag & INVOKED_SCALAR)) {
      temperature->compute_scalar();
      temperature->invoked_flag |= INVOKED_SCALAR;
    }
    compute_ke();

  } else if (strcmp(word,"etotal") == 0) {
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    if (update->whichflag == 0) {
      if (pe->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pe->invoked_flag & INVOKED_SCALAR)) {
      pe->compute_scalar();
      pe->invoked_flag |= INVOKED_SCALAR;
    }
    if (!temperature)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init temp");
    if (update->whichflag == 0) {
      if (temperature->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(temperature->invoked_flag & INVOKED_SCALAR)) {
      temperature->compute_scalar();
      temperature->invoked_flag |= INVOKED_SCALAR;
    }
    compute_etotal();

  } else if (strcmp(word,"enthalpy") == 0) {
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    if (update->whichflag == 0) {
      if (pe->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pe->invoked_flag & INVOKED_SCALAR)) {
      pe->compute_scalar();
      pe->invoked_flag |= INVOKED_SCALAR;
    }
    if (!temperature)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init temp");
    if (update->whichflag == 0) {
      if (temperature->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(temperature->invoked_flag & INVOKED_SCALAR)) {
      temperature->compute_scalar();
      temperature->invoked_flag |= INVOKED_SCALAR;
    }
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_scalar != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_SCALAR)) {
      pressure->compute_scalar();
      pressure->invoked_flag |= INVOKED_SCALAR;
    }
    compute_enthalpy();

  } else if (strcmp(word,"evdwl") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_evdwl();

  } else if (strcmp(word,"ecoul") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_ecoul();

  } else if (strcmp(word,"epair") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_epair();

  } else if (strcmp(word,"ebond") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_ebond();

  } else if (strcmp(word,"eangle") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_eangle();

  } else if (strcmp(word,"edihed") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_edihed();

  } else if (strcmp(word,"eimp") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_eimp();

  } else if (strcmp(word,"emol") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_emol();

  } else if (strcmp(word,"elong") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_elong();

  } else if (strcmp(word,"etail") == 0) {
    if (update->eflag_global != update->ntimestep)
      error->all("Energy was not tallied on needed timestep");
    if (!pe)
      error->all("Thermo keyword in variable requires thermo to use/init pe");
    pe->invoked_flag |= INVOKED_SCALAR;
    compute_etail();

  } else if (strcmp(word,"vol") == 0) compute_vol();
  else if (strcmp(word,"lx") == 0) compute_lx();
  else if (strcmp(word,"ly") == 0) compute_ly();
  else if (strcmp(word,"lz") == 0) compute_lz();

  else if (strcmp(word,"xlo") == 0) compute_xlo();
  else if (strcmp(word,"xhi") == 0) compute_xhi();
  else if (strcmp(word,"ylo") == 0) compute_ylo();
  else if (strcmp(word,"yhi") == 0) compute_yhi();
  else if (strcmp(word,"zlo") == 0) compute_zlo();
  else if (strcmp(word,"zhi") == 0) compute_zhi();

  else if (strcmp(word,"xy") == 0) compute_xy();
  else if (strcmp(word,"xz") == 0) compute_xz();
  else if (strcmp(word,"yz") == 0) compute_yz();

  else if (strcmp(word,"pxx") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pxx();

  } else if (strcmp(word,"pyy") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pyy();

  } else if (strcmp(word,"pzz") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pzz();

  } else if (strcmp(word,"pxy") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pxy();

  } else if (strcmp(word,"pxz") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pxz();

  } else if (strcmp(word,"pyz") == 0) {
    if (!pressure)
      error->all("Thermo keyword in variable requires "
		 "thermo to use/init press");
    if (update->whichflag == 0) {
      if (pressure->invoked_vector != update->ntimestep)
	error->all("Compute used in variable thermo keyword between runs "
		   "is not current");
    } else if (!(pressure->invoked_flag & INVOKED_VECTOR)) {
      pressure->compute_vector();
      pressure->invoked_flag |= INVOKED_VECTOR;
    }
    compute_pyz();

  } else return 1;

  *answer = dvalue;
  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
   compute/fix are normalized by atoms if returning extensive value
   variable value is not normalized (formula should normalize if desired)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Thermo::compute_compute()
{
  Compute *compute = computes[field2index[ifield]];
  if (argindex[ifield] == 0) {
    dvalue = compute->scalar;
    if (normflag && compute->extscalar) dvalue /= natoms;
  } else {
    dvalue = compute->vector[argindex[ifield]-1];
    if (normflag) {
      if (compute->extvector == 0) return;
      else if (compute->extvector == 1) dvalue /= natoms;
      else if (compute->extlist[argindex[ifield]-1]) dvalue /= natoms;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fix()
{
  Fix *fix = fixes[field2index[ifield]];
  if (argindex[ifield] == 0) {
    dvalue = fix->compute_scalar();
    if (normflag && fix->extscalar) dvalue /= natoms;
  } else {
    dvalue = fix->compute_vector(argindex[ifield]-1);
    if (normflag) {
      if (fix->extvector == 0) return;
      else if (fix->extvector == 1) dvalue /= natoms;
      else if (fix->extlist[argindex[ifield]-1]) dvalue /= natoms;
    }
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_variable()
{
  dvalue = input->variable->compute_equal(variables[field2index[ifield]]);
}

/* ----------------------------------------------------------------------
   one method for every keyword thermo can output
   called by compute() or evaluate_keyword()
   compute will have already been called
   set ivalue/dvalue if value is integer/double
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void Thermo::compute_step()
{
  ivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_atoms()
{
  ivalue = static_cast<int> (natoms);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cpu()
{
  if (firststep == 0) dvalue = 0.0;
  else dvalue = timer->elapsed(TIME_LOOP);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_temp()
{
  dvalue = temperature->scalar;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_press()
{
  dvalue = pressure->scalar;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pe()
{
  dvalue = pe->scalar;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ke()
{
  dvalue = temperature->scalar;
  dvalue *= 0.5 * temperature->dof * force->boltz;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etotal()
{
  compute_pe();
  double ke = temperature->scalar;
  ke *= 0.5 * temperature->dof * force->boltz;
  if (normflag) ke /= natoms;
  dvalue += ke;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_enthalpy()
{
  compute_etotal();
  double etmp = dvalue;

  compute_vol();
  double vtmp = dvalue;
  if (normflag) vtmp /= natoms;

  compute_press();
  double ptmp = dvalue;

  dvalue = etmp + ptmp*vtmp/(force->nktv2p);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_evdwl()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ecoul()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_coul;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_epair()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) dvalue += force->kspace->energy;
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ebond()
{
  if (force->bond) {
    double tmp = force->bond->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eangle()
{
  if (force->angle) {
    double tmp = force->angle->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_edihed()
{
  if (force->dihedral) {
    double tmp = force->dihedral->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eimp()
{
  if (force->improper) {
    double tmp = force->improper->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_emol()
{
  double tmp = 0.0;
  if (atom->molecular) {
    if (force->bond) tmp += force->bond->energy;
    if (force->angle) tmp += force->angle->energy;
    if (force->dihedral) tmp += force->dihedral->energy;
    if (force->improper) tmp += force->improper->energy;
    MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_elong()
{
  if (force->kspace) {
    dvalue = force->kspace->energy;
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etail()
{
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue = force->pair->etail / volume;
    if (normflag) dvalue /= natoms;
  } else dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_vol()
{
  if (domain->dimension == 3)
    dvalue = domain->xprd * domain->yprd * domain->zprd;
  else
    dvalue = domain->xprd * domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lx()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ly()
{
  dvalue = domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lz()
{
  dvalue = domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xlo()
{
  dvalue = domain->boxlo[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xhi()
{
  dvalue = domain->boxhi[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ylo()
{
  dvalue = domain->boxlo[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_yhi()
{
  dvalue = domain->boxhi[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zlo()
{
  dvalue = domain->boxlo[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zhi()
{
  dvalue = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xy()
{
  dvalue = domain->xy;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xz()
{
  dvalue = domain->xz;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_yz()
{
  dvalue = domain->yz;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxx()
{
  dvalue = pressure->vector[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyy()
{
  dvalue = pressure->vector[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pzz()
{
  dvalue = pressure->vector[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxy()
{
  dvalue = pressure->vector[3];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxz()
{
  dvalue = pressure->vector[4];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyz()
{
  dvalue = pressure->vector[5];
}
