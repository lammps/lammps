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
#include "string.h"
#include "stdlib.h"
#include "dump_local.h"
#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{INT,DOUBLE};

#define INVOKED_LOCAL 16

/* ---------------------------------------------------------------------- */

DumpLocal::DumpLocal(LAMMPS *lmp, int narg, char **arg) : 
  Dump(lmp, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump local arguments specified");

  clearstep = 1;

  nevery = atoi(arg[3]);

  size_one = nfield = narg-5;
  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  // computes & fixes which the dump accesses

  field2index = new int[nfield];
  argindex = new int[nfield];

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  // process attributes

  parse_fields(narg,arg);

  // setup format strings

  vformat = new char*[size_one];

  format_default = new char[3*size_one+1];
  format_default[0] = '\0';

  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == INT) format_default = strcat(format_default,"%d ");
    else format_default = strcat(format_default,"%g ");
    vformat[i] = NULL;
  }

  // setup column string

  int n = 0;
  for (int iarg = 5; iarg < narg; iarg++) n += strlen(arg[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 5; iarg < narg; iarg++) {
    strcat(columns,arg[iarg]);
    strcat(columns," ");
  }

  // setup default label string

  char *str = (char *) "ENTRIES";
  n = strlen(str) + 1;
  label = new char[n];
  strcpy(label,str);
}

/* ---------------------------------------------------------------------- */

DumpLocal::~DumpLocal()
{
  delete [] pack_choice;
  delete [] vtype;
  delete [] field2index;
  delete [] argindex;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
  delete [] label;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::init_style()
{
  if (sort_flag && sortcol == 0)
    error->all(FLERR,"Dump local cannot sort by atom ID");

  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    if (i == 0) ptr = strtok(format," \0");
    else ptr = strtok(NULL," \0");
    delete [] vformat[i];
    vformat[i] = new char[strlen(ptr) + 2];
    strcpy(vformat[i],ptr);
    vformat[i] = strcat(vformat[i]," ");
  }

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump local compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump local fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->local_freq)
      error->all(FLERR,"Dump local and fix not computed at compatible times");
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpLocal::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"label") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    delete [] label;
    int n = strlen(arg[1]) + 1;
    label = new char[n];
    strcpy(label,arg[1]);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_header(bigint ndump)
{
  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF %s\n",label);
    fprintf(fp,BIGINT_FORMAT "\n",ndump);
    fprintf(fp,"ITEM: %s %s\n",label,columns);
  }
}

/* ---------------------------------------------------------------------- */

int DumpLocal::count()
{
  int i;

  // invoke Computes for local quantities

  if (ncompute) {
    for (i = 0; i < ncompute; i++) {
      if (!(compute[i]->invoked_flag & INVOKED_LOCAL)) {
	compute[i]->compute_local();
	compute[i]->invoked_flag |= INVOKED_LOCAL;
      }
    }
  }

  // nmine = # of local values I contribute
  // must be consistent for all input fields

  nmine = -1;

  for (int i = 0; i < ncompute; i++) {
    if (nmine < 0) nmine = compute[i]->size_local_rows;
    else if (nmine != compute[i]->size_local_rows)
      error->one(FLERR,"Dump local count is not consistent across input fields");
  }

  for (int i = 0; i < nfix; i++) {
    if (nmine < 0) nmine = fix[i]->size_local_rows;
    else if (nmine != fix[i]->size_local_rows)
      error->one(FLERR,"Dump local count is not consistent across input fields");
  }

  return nmine;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::pack(int *dummy)
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_data(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else fprintf(fp,vformat[j],mybuf[m]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocal::parse_fields(int narg, char **arg)
{
  int computefixflag = 0;

  // customize by adding to if statement

  int i;
  for (int iarg = 5; iarg < narg; iarg++) {
    i = iarg-5;

    if (strcmp(arg[iarg],"index") == 0) {
      pack_choice[i] = &DumpLocal::pack_index;
      vtype[i] = INT;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      computefixflag = 1;
      pack_choice[i] = &DumpLocal::pack_compute;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump local command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump local compute ID");
      if (modify->compute[n]->local_flag == 0)
	error->all(FLERR,"Dump local compute does not compute local info");
      if (argindex[i] == 0 && modify->compute[n]->size_local_cols > 0)
	error->all(FLERR,"Dump local compute does not calculate local vector");
      if (argindex[i] > 0 && modify->compute[n]->size_local_cols == 0)
	error->all(FLERR,"Dump local compute does not calculate local array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->compute[n]->size_local_cols)
	error->all(FLERR,"Dump local compute vector is accessed out-of-range");

      field2index[i] = add_compute(suffix);
      delete [] suffix;
      
    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      computefixflag = 1;
      pack_choice[i] = &DumpLocal::pack_fix;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Invalid attribute in dump local command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump local fix ID");
      if (modify->fix[n]->local_flag == 0)
	error->all(FLERR,"Dump local fix does not compute local info");
      if (argindex[i] == 0 && modify->fix[n]->size_local_cols > 0)
	error->all(FLERR,"Dump local fix does not compute local vector");
      if (argindex[i] > 0 && modify->fix[n]->size_local_cols == 0)
	error->all(FLERR,"Dump local fix does not compute local array");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->fix[n]->size_local_cols)
	error->all(FLERR,"Dump local fix vector is accessed out-of-range");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    } else error->all(FLERR,"Invalid attribute in dump local command");
  }

  if (computefixflag == 0)
    error->all(FLERR,"Dump local attributes contain no compute or fix");
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpLocal::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;
  
  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpLocal::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;
  
  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix results
------------------------------------------------------------------------- */

void DumpLocal::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->vector_local;
  double **array = compute[field2index[n]]->array_local;
  int ncount = compute[field2index[n]]->size_local_rows;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < ncount; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < ncount; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpLocal::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->vector_local;
  double **array = fix[field2index[n]]->array_local;
  int index = argindex[n];

  if (index == 0) {
    for (int i = 0; i < nmine; i++) {
      buf[n] = vector[i];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nmine; i++) {
      buf[n] = array[i][index];
      n += size_one;
    }
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump local can output
   the local value is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void DumpLocal::pack_index(int n)
{
  int index;
  MPI_Scan(&nmine,&index,1,MPI_INT,MPI_SUM,world);
  index -= nmine;

  for (int i = 0; i < nmine; i++) {
    buf[n] = ++index;
    n += size_one;
  }
}
