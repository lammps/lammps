/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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
#include "force.h"
#include "temperature.h"
#include "pressure.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
//#include "pair_dipole_cut.h"
//#include "pair_dipole_long.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "fix.h"
#include "comm.h"
#include "atom.h"
#include "domain.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

// customize by adding a keyword to this list:

// step, atoms, cpu, temp, press, pe, ke, eng,
// evdwl, ecoul, epair, ebond, eangle, edihed, eimp, emol, elong, etail, erot
// vol, lx, ly, lz, pxx, pyy, pzz, pxy, pxz, pyz
// gke, grot (granular trans-ke and rotational-ke)
// tave, pave, eave, peave (time-averaged quantities)
// t_ID (user-defined temperatures)

// customize by adding a DEFINE to this list

#define ONE "step temp epair emol eng press"
#define MULTI "eng ke temp pe ebond eangle edihed eimp evdwl ecoul elong press"
#define GRANULAR "step atoms gke grot"
#define DIPOLE "step temp epair erot eng press"

enum {IGNORE,WARN,ERROR};           // same as write_restart.cpp
enum {ONELINE,MULTILINE};
enum {INT,FLOAT};
enum {NEITHER,PRINT,ENERGY,BOTH};   // same as modify.cpp

#define MAXLINE 1024
#define INERTIA3D 0.4               // moments of inertia for sphere and disk
#define INERTIA2D 0.5

/* ---------------------------------------------------------------------- */

Thermo::Thermo(int narg, char **arg)
{
  MPI_Comm_rank(world,&me);

  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  // set thermo_modify defaults

  normuserflag = 0;
  lineflag = ONELINE;
  lostflag = ERROR;
  lostbefore = 0;
  flushflag = 0;
  nwindow = 10;
  ncurrent_t = ncurrent_p = ncurrent_e = ncurrent_pe = -1;
  npartial_t = npartial_p = npartial_e = npartial_pe = 0;

  // initialize ptrs

  // set style and corresponding lineflag
  // custom style builds its own line of keywords
  // customize a new DEFINE style by adding to if statement

  line = new char[MAXLINE];

  if (strcmp(style,"one") == 0) {
    strcpy(line,ONE);
  } else if (strcmp(style,"multi") == 0) {
    strcpy(line,MULTI);
    lineflag = MULTILINE;
  } else if (strcmp(style,"granular") == 0) {
    strcpy(line,GRANULAR);
  } else if (strcmp(style,"dipole") == 0) {
    strcpy(line,DIPOLE);
  } else if (strcmp(style,"custom") == 0) {
    if (narg == 1) error->all("Illegal thermo style custom command");
    line[0] = '\0';
    for (int iarg = 1; iarg < narg; iarg++) {
      strcat(line,arg[iarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';
  } else error->all("Illegal thermo style command");

  // process line of keywords

  nfield = ntemp = nfix = 0;
  tempindices = NULL;
  parse_fields(line);
  nfield_initial = nfield;

  // format strings

  char *str = "%8d";
  n = strlen(str) + 1;
  format_int_def = new char[n];
  strcpy(format_int_def,str);

  str = "%12.8g";
  n = strlen(str) + 1;
  format_g_def = new char[n];
  strcpy(format_g_def,str);

  str = "%14.4f";
  n = strlen(str) + 1;
  format_f_def = new char[n];
  strcpy(format_f_def,str);

  str = "---------------- Step %8d ----- CPU = %11.4f (sec) ----------------";
  n = strlen(str) + 1;
  format_multi = new char[n];
  strcpy(format_multi,str);

  format_iflag = 0;
  format_fflag = 0;
  format_int_user = NULL;
  format_float_user = NULL;

  nformat = 0;
  format_user = NULL;
  format_index = NULL;

  format = NULL;

  // ptrs for T,P computation
  // temperature is default in force->templist

  temperature = force->templist[0];
  pressure = force->pressure;

  // average quantities

  tsum = psum = esum = pesum = 0.0;
  tpast = new double[nwindow];
  ppast = new double[nwindow];
  epast = new double[nwindow];
  pepast = new double[nwindow];

  inertia = NULL;
}

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
  delete [] style;
  delete [] line;

  if (keyword) {
    for (int i = 0; i < nfield; i++) delete [] keyword[i];
    memory->sfree(keyword);
    memory->sfree(vfunc);
    memory->sfree(vtype);
  }

  delete [] format_int_def;
  delete [] format_g_def;
  delete [] format_f_def;
  delete [] format_multi;
  if (format_iflag) delete [] format_int_user;
  if (format_fflag) delete [] format_float_user;

  if (nformat) {
    for (int i = 0; i < nformat; i++) delete [] format_user[i];
    memory->sfree(format_user);
    memory->sfree(format_index);
  }

  if (format) {
    for (int i = 0; i < nfield; i++) delete [] format[i];
    delete [] format;
  }

  if (tempindices) memory->sfree(tempindices);

  if (nfix) {
    delete [] fixflags;
    delete [] fixprint;
    delete [] fixenergy;
    delete [] fixvalues;
  }

  delete [] tpast;
  delete [] ppast;
  delete [] epast;
  delete [] pepast;

  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
  int i,j,n;

  // error checks
  // mass/rmass for granular and non-granular styles
  // dipole/omega for erot

  if (granflag && atom->check_style("granular") == 0)
    error->all("Must use atom style granular with chosen thermo settings");
  if ((tempflag || pressflag) && atom->check_style("granular"))
    error->all("Cannot use atom style granular with chosen thermo settings");
  if (dipoleflag && atom->check_style("dipole") == 0)
    error->all("Must use atom style dipole with chosen thermo settings");

  // set normvalue to default setting unless user has specified it

  if (normuserflag) normvalue = normuser;
  else if (strcmp(style,"granular") == 0) normvalue = 0;
  else if (strcmp(update->unit_style,"lj") == 0) normvalue = 1;
  else normvalue = 0;

  // check for granular freeze Fix and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // setup moments of inertia if dipoleflag is set
  // extract sigma as particle size to compute moment of inertia

  if (dipoleflag) {
    delete [] inertia;
    inertia = new double[atom->ntypes+1];
    double *mass = atom->mass;

    Pair *anypair;
    double **sigma;
    /*
    if (anypair = force->pair_match("dipole/cut"))
      sigma = ((PairDipoleCut *) anypair)->sigma;
    else if (anypair = force->pair_match("dipole/long"))
      sigma = ((PairDipoleLong *) anypair)->sigma;
    else error->all("Pair style is incompatible with thermo_style erot");
    */

    if (force->dimension == 3)
      for (int i = 1; i <= atom->ntypes; i++)
	inertia[i] = INERTIA3D * mass[i] * 0.25*sigma[i][i]*sigma[i][i];
    else
      for (int i = 1; i <= atom->ntypes; i++)
	inertia[i] = INERTIA2D * mass[i] * 0.25*sigma[i][i]*sigma[i][i];
  }

  // clear all format strings
  // clear extra keyword strings beyond nfield_initial

  if (format) {
    for (i = 0; i < nfield; i++) delete [] format[i];
    delete [] format;
  }
  for (i = nfield_initial; i < nfield; i++) delete [] keyword[i];
  nfield = nfield_initial;

  // add a field if volume changes and not style = custom
  // this check must come after domain init, so box_change is set

  if (domain->box_change && strcmp(style,"custom") != 0)
    addfield("Volume",&Thermo::compute_vol,FLOAT);
  
  // add fields due to fixes
  
  if (nfix) {
    delete [] fixflags;
    delete [] fixprint;
    delete [] fixenergy;
    delete [] fixvalues;
  }

  nfix = nfix_energy = nfix_print = 0;
  if (modify->n_thermo) nfix = modify->thermo_fields(0,NULL,NULL);
  if (nfix) {
    fixflags = new int[nfix];
    fixprint = new int[nfix];
    fixenergy = new int[nfix];
    fixvalues = new double[nfix];
    char **strings = new char*[nfix];
    for (i = 0; i < nfix; i++) strings[i] = new char[9];

    nfix = modify->thermo_fields(nfix,fixflags,strings);

    for (i = 0; i < nfix; i++) {
      if (fixflags[i] == PRINT || fixflags[i] == BOTH) {
	fixprint[nfix_print++] = i;
	addfield(strings[i],&Thermo::compute_fix,FLOAT);
      }
      if (fixflags[i] == ENERGY || fixflags[i] == BOTH)
	fixenergy[nfix_energy++] = i;
    }

    for (i = 0; i < nfix; i++) delete [] strings[i];
    delete [] strings;
  }

  // nextra = 1 if printing extra temperature or fix values

  nextra = 0;
  if (ntemp || nfix) nextra = 1;

  // set format string for each printed value
  // include keyword if lineflag = MULTILINE
  // add '/n' every 3 values if lineflag = MULTILINE
  // add trailing '/n' to last value

  format = new char*[nfield];
  for (int i = 0; i < nfield; i++) format[i] = new char[32];
  char *ptr;

  for (i = 0; i < nfield; i++) {
    format[i][0] = '\0';
    if (lineflag == MULTILINE && i % 3 == 0) strcat(format[i],"\n");

    for (j = nformat-1; j >= 0; j--)
      if (format_index[j] == i+1) break;
    if (j >= 0) ptr = format_user[j];
    else if (vtype[i] == INT && format_iflag) ptr = format_int_user;
    else if (vtype[i] == INT) ptr = format_int_def;
    else if (vtype[i] == FLOAT && format_fflag) ptr = format_float_user;
    else if (lineflag == MULTILINE) ptr = format_f_def;
    else ptr = format_g_def;

    n = strlen(format[i]);
    if (lineflag == ONELINE) sprintf(&format[i][n],"%s ",ptr);
    else sprintf(&format[i][n],"%-8s = %s ",keyword[i],ptr);

    if (i == nfield-1) strcat(format[i],"\n");
  }
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
  firststep = flag;

  // check for lost atoms
  // turn off normflag if natoms = 0 to avoid divide by 0
  // compute requested thermo quantities

  natoms = lost_check();
  if (natoms == 0) normflag = 0;
  else normflag = normvalue;

  if (tempflag) tempvalue = temperature->compute();
  if (pressflag) pressure->compute(temperature);
  if (nextra) {
    itemp_print = ifix_print = 0;
    if (nfix) modify->thermo_compute(fixvalues);
  }

  // if lineflag = MULTILINE, prepend step/cpu header line

  int loc = 0;
  if (lineflag == MULTILINE) {
    double cpu;
    if (flag) cpu = timer->elapsed(TIME_LOOP);
    else cpu = 0.0;
    loc = sprintf(&line[loc],format_multi,update->ntimestep,cpu);
  }

  // add each thermo value to line with its specific format

  for (int i = 0; i < nfield; i++) {
    (this->*vfunc[i])();
    if (vtype[i] == INT) loc += sprintf(&line[loc],format[i],ivalue);
    else loc += sprintf(&line[loc],format[i],dvalue);
  }

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
  // natoms = current # of atoms

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);
  if (natoms == atom->natoms) return natoms;

  // if not checking or already warned, just return

  if (lostflag == IGNORE) return natoms;
  if (lostflag == WARN && lostbefore == 1) return natoms;

  // error message

  if (lostflag == ERROR) {
    char str[128];
    sprintf(str,"Lost atoms: original %.15g current %.15g",
	    atom->natoms,natoms);
    error->all(str);
  }

  // warning message

  char str[128];
  sprintf(str,"Lost atoms: original %.15g current %.15g",atom->natoms,natoms);
  if (me == 0) error->warning(str);
  lostbefore = 1;
  return natoms;
}

/* ----------------------------------------------------------------------
   modify thermo parameters
------------------------------------------------------------------------- */

void Thermo::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal thermo_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      int tempwhich;
      for (tempwhich = 0; tempwhich < force->ntemp; tempwhich++)
	if (strcmp(arg[iarg+1],force->templist[tempwhich]->id) == 0) break;
      if (tempwhich == force->ntemp) 
	error->all("Could not find thermo_modify temperature ID");
      if (force->templist[tempwhich]->igroup != 0 && comm->me == 0)
	error->warning("Temperature for thermo pressure is not for group all");
      if (strcmp(force->templist[tempwhich]->style,"region") == 0 && 
	  comm->me == 0)
	error->warning("Temperature for thermo pressure is style region");
      temperature = force->templist[tempwhich];
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
	if (format_iflag) delete [] format_int_user;
	format_iflag = 1;
	int n = strlen(arg[iarg+2]) + 1;
	format_int_user = new char[n];
	strcpy(format_int_user,arg[iarg+2]);
      } else if (strcmp(arg[iarg+1],"float") == 0) {
	if (format_fflag) delete [] format_float_user;
	format_fflag = 1;
	int n = strlen(arg[iarg+2]) + 1;
	format_float_user = new char[n];
	strcpy(format_float_user,arg[iarg+2]);
      } else {
	int i = atoi(arg[iarg+1]);
	if (i < 1) error->all("Illegal thermo_modify command");
	format_user =
	  (char **) memory->srealloc(format_user,(nformat+1)*sizeof(char *),
				     "thermo:format_user");
	format_index =
	  (int *) memory->srealloc(format_index,(nformat+1)*sizeof(int),
				   "thermo:format_index");
	int n = strlen(arg[iarg+2]) + 1;
	format_user[nformat] = new char[n];
	strcpy(format_user[nformat],arg[iarg+2]);
	format_index[nformat] = i;
	nformat++;
      }
      iarg += 3;
    } else if (strcmp(arg[iarg],"window") == 0) {
      if (iarg+2 > narg) error->all("Illegal thermo_modify command");
      nwindow = atoi(arg[iarg+1]);
      if (nwindow <= 0) error->all("Illegal thermo_modify command");
      ncurrent_t = ncurrent_p = ncurrent_e = ncurrent_pe = -1;
      npartial_t = npartial_p = npartial_e = npartial_pe = 0;
      tsum = psum  = esum = pesum = 0.0;
      delete [] tpast;
      delete [] ppast;
      delete [] epast;
      delete [] pepast;
      tpast = new double[nwindow];
      ppast = new double[nwindow];
      epast = new double[nwindow];
      pepast = new double[nwindow];
      iarg += 2;
    } else error->all("Illegal thermo_modify command");
  }
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(char *str)
{
  keyword = NULL;
  vfunc = NULL;
  vtype = NULL;

  tempflag = pressflag = tensorflag = granflag = peflag = dipoleflag = 0;

  // customize by adding to if statement

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
      tempflag = 1;
    } else if (strcmp(word,"press") == 0) {
      addfield("Press",&Thermo::compute_press,FLOAT);
      tempflag = pressflag = 1;
    } else if (strcmp(word,"pe") == 0) {
      peflag = 1;
      addfield("PotEng",&Thermo::compute_pe,FLOAT);
    } else if (strcmp(word,"ke") == 0) {
      addfield("KinEng",&Thermo::compute_ke,FLOAT);
      tempflag = 1;
    } else if (strcmp(word,"eng") == 0) {
      addfield("TotEng",&Thermo::compute_eng,FLOAT);
      tempflag = 1;
    } else if (strcmp(word,"evdwl") == 0) {
      addfield("E_vdwl",&Thermo::compute_evdwl,FLOAT);
    } else if (strcmp(word,"ecoul") == 0) {
      addfield("E_coul",&Thermo::compute_ecoul,FLOAT);
    } else if (strcmp(word,"epair") == 0) {
      addfield("E_pair",&Thermo::compute_epair,FLOAT);
    } else if (strcmp(word,"ebond") == 0) {
      addfield("E_bond",&Thermo::compute_ebond,FLOAT);
    } else if (strcmp(word,"eangle") == 0) {
      addfield("E_angle",&Thermo::compute_eangle,FLOAT);
    } else if (strcmp(word,"edihed") == 0) {
      addfield("E_dihed",&Thermo::compute_edihed,FLOAT);
    } else if (strcmp(word,"eimp") == 0) {
      addfield("E_impro",&Thermo::compute_eimp,FLOAT);
    } else if (strcmp(word,"emol") == 0) {
      addfield("E_mol",&Thermo::compute_emol,FLOAT);
    } else if (strcmp(word,"elong") == 0) {
      addfield("E_long",&Thermo::compute_elong,FLOAT);
    } else if (strcmp(word,"etail") == 0) {
      addfield("E_tail",&Thermo::compute_etail,FLOAT);
    } else if (strcmp(word,"erot") == 0) {
      addfield("E_rot",&Thermo::compute_erot,FLOAT);
      dipoleflag = 1;
    } else if (strcmp(word,"vol") == 0) {
      addfield("Volume",&Thermo::compute_vol,FLOAT);
    } else if (strcmp(word,"lx") == 0) {
      addfield("Lx",&Thermo::compute_lx,FLOAT);
    } else if (strcmp(word,"ly") == 0) {
      addfield("Ly",&Thermo::compute_ly,FLOAT);
    } else if (strcmp(word,"lz") == 0) {
      addfield("Lz",&Thermo::compute_lz,FLOAT);
    } else if (strcmp(word,"pxx") == 0) {
      addfield("Pxx",&Thermo::compute_pxx,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"pyy") == 0) {
      addfield("Pyy",&Thermo::compute_pyy,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"pzz") == 0) {
      addfield("Pzz",&Thermo::compute_pzz,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"pxy") == 0) {
      addfield("Pxy",&Thermo::compute_pxy,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"pxz") == 0) {
      addfield("Pxz",&Thermo::compute_pxz,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"pyz") == 0) {
      addfield("Pyz",&Thermo::compute_pyz,FLOAT);
      tempflag = pressflag = tensorflag = 1;
    } else if (strcmp(word,"gke") == 0) {
      addfield("Trans-KE",&Thermo::compute_gke,FLOAT);
      granflag = 1;
    } else if (strcmp(word,"grot") == 0) {
      addfield("Rot-KE",&Thermo::compute_grot,FLOAT);
      granflag = 1;
    } else if (strcmp(word,"tave") == 0) {
      addfield("T_ave",&Thermo::compute_tave,FLOAT);
      tempflag = 1;
    } else if (strcmp(word,"pave") == 0) {
      addfield("P_ave",&Thermo::compute_pave,FLOAT);
      pressflag = 1;
    } else if (strcmp(word,"eave") == 0) {
      addfield("E_ave",&Thermo::compute_eave,FLOAT);
    } else if (strcmp(word,"peave") == 0) {
      addfield("PE_ave",&Thermo::compute_peave,FLOAT);

    // user-defined temperature (t_ID)
    // only 1st 8 chars of ID are passed to addfield
    // extend tempindices and store index of ID into force::templist

    } else if (strncmp(word,"t_",2) == 0) {
      char copy[9];
      strncpy(copy,word,8);
      copy[0] = 'T';
      copy[8] = '\0';
      addfield(copy,&Thermo::compute_t_id,FLOAT);
      int tempwhich;
      for (tempwhich = 0; tempwhich < force->ntemp; tempwhich++)
	if (strcmp(&word[2],force->templist[tempwhich]->id) == 0) break;
      if (tempwhich == force->ntemp) 
	error->all("Could not find thermo temperature ID");
      tempindices = (int *) memory->srealloc(tempindices,
					     (ntemp+1)*sizeof(int),
					     "thermo:tempindices");
      tempindices[ntemp++] = tempwhich;
    } else error->all("Invalid keyword in thermo_style command");

    word = strtok(NULL," \0");
  }
}

/* ----------------------------------------------------------------------
   add one field to list of quantities to print
   extend keyword,vfunc,vtype arrays
------------------------------------------------------------------------- */

void Thermo::addfield(char *key, FnPtr func, int typeflag)
{
  keyword = (char **)
    memory->srealloc(keyword,(nfield+1)*sizeof(char *),"thermo:keyword");
  vfunc = (FnPtr *)
    memory->srealloc(vfunc,(nfield+1)*sizeof(FnPtr),"thermo:vfunc");
  vtype = (int *)
    memory->srealloc(vtype,(nfield+1)*sizeof(int),"thermo:vtype");

  int n = strlen(key) + 1;
  keyword[nfield] = new char[n];
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   compute a single thermodyanmic value matching word to custom list
   called by variable evaulation function from input script
   return it as double in answer
   return 0 if OK, 1 if str is invalid
   customize by adding to if statement
     not all keywords are supported
     only add keyword if it makes sense in input script variable
------------------------------------------------------------------------- */

int Thermo::compute_value(char *word, double *answer)
{
  if (strcmp(word,"temp") == 0) {
    tempvalue = temperature->compute();
    compute_temp();
  } else if (strcmp(word,"press") == 0) {
    tempvalue = temperature->compute();
    pressure->compute(temperature);
    compute_press();
  } else if (strcmp(word,"pe") == 0) {
    fix_compute_pe();
    compute_pe();
  } else if (strcmp(word,"ke") == 0) {
    tempvalue = temperature->compute();
    compute_ke();
  } else if (strcmp(word,"eng") == 0) {
    tempvalue = temperature->compute();
    fix_compute_pe();
    compute_eng();
  } else if (strcmp(word,"evdwl") == 0) compute_evdwl();
  else if (strcmp(word,"ecoul") == 0) compute_ecoul();
  else if (strcmp(word,"epair") == 0) compute_epair();
  else if (strcmp(word,"ebond") == 0) compute_ebond();
  else if (strcmp(word,"eangle") == 0) compute_eangle();
  else if (strcmp(word,"edihed") == 0) compute_edihed();
  else if (strcmp(word,"eimp") == 0) compute_eimp();
  else if (strcmp(word,"emol") == 0) compute_emol();
  else if (strcmp(word,"elong") == 0) compute_elong();
  else if (strcmp(word,"etail") == 0) compute_etail();
  else if (strcmp(word,"erot") == 0) compute_erot();
  else if (strcmp(word,"gke") == 0) compute_gke();
  else if (strcmp(word,"grot") == 0) compute_grot();
  else if (strncmp(word,"t_",2) == 0) {
    int tempwhich;
    for (tempwhich = 0; tempwhich < force->ntemp; tempwhich++)
      if (strcmp(&word[2],force->templist[tempwhich]->id) == 0) break;
    if (tempwhich == force->ntemp) return 1;
    dvalue = force->templist[tempwhich]->compute();
  } else return 1;
  
  *answer = dvalue;
  return 0;
}

/* ----------------------------------------------------------------------
   compute potential energy from any fixes
   called when compute_pe() is called directly and not via Thermo::compute()
   e.g. from minimize or Thermo::compute_value()
------------------------------------------------------------------------- */

void Thermo::fix_compute_pe()
{
  if (nfix_energy) modify->thermo_compute(fixvalues);
}

/* ----------------------------------------------------------------------
   one routine for every quantity thermo can output
   set ivalue/dvalue if value is integer/double
   customize by adding a method
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
  dvalue = tempvalue;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_press()
{
  dvalue = pressure->p_total;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pe()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) 
    tmp += force->dihedral->eng_vdwl + force->dihedral->eng_coul;

  if (atom->molecular) {
    if (force->bond) tmp += force->bond->energy;
    if (force->angle) tmp += force->angle->energy;
    if (force->dihedral) tmp += force->dihedral->energy;
    if (force->improper) tmp += force->improper->energy;
  }

  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) dvalue += force->kspace->energy;
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }
  if (nfix_energy)
    for (int i = 0; i < nfix_energy; i++) dvalue += fixvalues[fixenergy[i]];

  if (normflag) dvalue /= natoms;

  // also store PE in potential_energy so other classes can grab it

  potential_energy = dvalue;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ke()
{
  dvalue = 0.5 * temperature->dof * force->boltz * tempvalue;
  if (normflag) dvalue /= natoms;
}


/* ---------------------------------------------------------------------- */

void Thermo::compute_eng()
{
  compute_pe();
  double ke = 0.5 * temperature->dof * force->boltz * tempvalue;
  if (normflag) ke /= natoms;
  dvalue += ke;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_evdwl()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) tmp += force->dihedral->eng_vdwl;

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
  if (force->dihedral) tmp += force->dihedral->eng_coul;
  MPI_Allreduce(&tmp,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_epair()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  if (force->bond) tmp += force->bond->eng_vdwl;
  if (force->dihedral) 
    tmp += force->dihedral->eng_vdwl + force->dihedral->eng_coul;
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

void Thermo::compute_erot()
{
  int *type = atom->type;
  double *dipole = atom->dipole;
  double **omega = atom->omega;
  int nlocal = atom->nlocal;

  double erot = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (dipole[type[i]] > 0.0)
      erot += inertia[type[i]] * 
	(omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + 
	 omega[i][2]*omega[i][2]);
  }
  double all;
  MPI_Allreduce(&erot,&all,1,MPI_DOUBLE,MPI_SUM,world);
  dvalue = 0.5 * all;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_vol()
{
  dvalue = domain->xprd * domain->yprd * domain->zprd;
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

void Thermo::compute_pxx()
{
  dvalue = pressure->p_tensor[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyy()
{
  dvalue = pressure->p_tensor[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pzz()
{
  dvalue = pressure->p_tensor[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxy()
{
  dvalue = pressure->p_tensor[3];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxz()
{
  dvalue = pressure->p_tensor[4];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyz()
{
  dvalue = pressure->p_tensor[5];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_gke()
{
  double *rmass = atom->rmass;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double nactive = 0.0;
  double ke = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & freeze_group_bit)) {
      nactive++;
      ke += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
    }

  ke *= 0.5 * force->mvv2e;
  MPI_Allreduce(&ke,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (normflag) {
    ke = dvalue;
    MPI_Allreduce(&nactive,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (dvalue > 0) dvalue = ke / dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_grot()
{
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double **phiv = atom->phiv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double nactive = 0.0;
  double ephi = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (!(mask[i] & freeze_group_bit)) {
      nactive++;
      ephi += (phiv[i][0]*phiv[i][0] + phiv[i][1]*phiv[i][1] + 
	       phiv[i][2]*phiv[i][2]) * radius[i]*radius[i]*rmass[i];
    }

  if (force->dimension == 3) ephi *= 0.5 * force->mvv2e * INERTIA3D;
  else ephi *= 0.5 * force->mvv2e * INERTIA2D;
  MPI_Allreduce(&ephi,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
  if (normflag) {
    ephi = dvalue;
    MPI_Allreduce(&nactive,&dvalue,1,MPI_DOUBLE,MPI_SUM,world);
    if (dvalue > 0) dvalue = ephi / dvalue;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_tave()
{
  if (firststep == 0) {
    if (npartial_t == 0) {
      ncurrent_t = 0;
      npartial_t = 1;
      tsum = tpast[0] = tempvalue;
    }
    dvalue = tsum/npartial_t;
  } else {
    ncurrent_t++;
    if (ncurrent_t == nwindow) ncurrent_t = 0;
    if (npartial_t == nwindow) tsum -= tpast[ncurrent_t];
    else npartial_t++;
    tpast[ncurrent_t] = tempvalue;
    tsum += tpast[ncurrent_t];
    dvalue = tsum/npartial_t;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pave()
{
  if (firststep == 0) {
    if (npartial_p == 0) {
      ncurrent_p = 0;
      npartial_p = 1;
      psum = ppast[0] = pressure->p_total;
    }
    dvalue = psum/npartial_p;
  } else {
    ncurrent_p++;
    if (ncurrent_p == nwindow) ncurrent_p = 0;
    if (npartial_p == nwindow) psum -= ppast[ncurrent_p];
    else npartial_p++;
    ppast[ncurrent_p] = pressure->p_total;
    psum += ppast[ncurrent_p];
    dvalue = psum/npartial_p;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eave()
{
  compute_eng();

  if (firststep == 0) {
    if (npartial_e == 0) {
      ncurrent_e = 0;
      npartial_e = 1;
      esum = epast[0] = dvalue;
    }
    dvalue = esum/npartial_e;
  } else {
    ncurrent_e++;
    if (ncurrent_e == nwindow) ncurrent_e = 0;
    if (npartial_e == nwindow) esum -= epast[ncurrent_e];
    else npartial_e++;
    epast[ncurrent_e] = dvalue;
    esum += epast[ncurrent_e];
    dvalue = esum/npartial_e;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_peave()
{
  compute_pe();

  if (firststep == 0) {
    if (npartial_pe == 0) {
      ncurrent_pe = 0;
      npartial_pe = 1;
      pesum = pepast[0] = dvalue;
    }
    dvalue = pesum/npartial_pe;
  } else {
    ncurrent_pe++;
    if (ncurrent_pe == nwindow) ncurrent_pe = 0;
    if (npartial_pe == nwindow) pesum -= pepast[ncurrent_pe];
    else npartial_pe++;
    pepast[ncurrent_pe] = dvalue;
    pesum += pepast[ncurrent_pe];
    dvalue = pesum/npartial_pe;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_t_id()
{
  dvalue = force->templist[tempindices[itemp_print++]]->compute();
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fix()
{
  dvalue = fixvalues[fixprint[ifix_print++]];
  if (normflag) dvalue /= natoms;
}
