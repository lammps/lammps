// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_local.h"

#include "arg_info.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpLocal::DumpLocal(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg),
  label(nullptr), vtype(nullptr), vformat(nullptr), columns(nullptr), field2index(nullptr),
  argindex(nullptr), id_compute(nullptr), compute(nullptr), id_fix(nullptr), fix(nullptr),
  pack_choice(nullptr)
{
  if (narg == 5) error->all(FLERR,"No dump local arguments specified");

  clearstep = 1;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal dump local command");

  if (binary)
    error->all(FLERR,"Binary files are not supported with dump local");

  nfield = narg - 5;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  nfield = utils::expand_args(FLERR,nfield,&arg[5],1,earg,lmp);

  if (earg != &arg[5]) expand = 1;

  // allocate field vectors

  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  buffer_allow = 1;
  buffer_flag = 1;

  // computes & fixes which the dump accesses

  field2index = new int[nfield];
  argindex = new int[nfield];

  ncompute = 0;
  id_compute = nullptr;
  compute = nullptr;

  nfix = 0;
  id_fix = nullptr;
  fix = nullptr;

  // process attributes

  parse_fields(nfield,earg);
  size_one = nfield;

  // setup format strings

  vformat = new char*[size_one];
  std::string fdefault;
  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == Dump::INT) fdefault += "%d ";
    else if (vtype[i] == Dump::DOUBLE) fdefault += "%g ";
    vformat[i] = nullptr;
  }
  format_default = utils::strdup(fdefault);

  format_column_user = new char*[size_one];
  for (int i = 0; i < size_one; i++) format_column_user[i] = nullptr;

  // setup column string

  std::string cols;
  for (int iarg = 0; iarg < nfield; iarg++) {
    cols += earg[iarg];
    cols += " ";
  }
  columns = utils::strdup(cols);

  // setup default label string

  label = utils::strdup("ENTRIES");

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nfield; i++) delete[] earg[i];
    memory->sfree(earg);
  }
}

/* ---------------------------------------------------------------------- */

DumpLocal::~DumpLocal()
{
  delete[] pack_choice;
  delete[] vtype;
  delete[] field2index;
  delete[] argindex;

  for (int i = 0; i < ncompute; i++) delete[] id_compute[i];
  memory->sfree(id_compute);
  delete[] compute;

  for (int i = 0; i < nfix; i++) delete[] id_fix[i];
  memory->sfree(id_fix);
  delete[] fix;

  for (int i = 0; i < size_one; i++) delete[] vformat[i];
  delete[] vformat;

  for (int i = 0; i < size_one; i++) delete[] format_column_user[i];
  delete[] format_column_user;

  delete[] columns;
  delete[] label;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::init_style()
{
  if (sort_flag && sortcol == 0)
    error->all(FLERR,"Dump local cannot sort by atom ID");

  // format = copy of default or user-specified line format

  delete[] format;
  if (format_line_user) format = utils::strdup(format_line_user);
  else format = utils::strdup(format_default);

  // tokenize the format string and add space at end of each format element
  // if user-specified int/float format exists, use it instead
  // if user-specified column format exists, use it instead
  // lo priority = line, medium priority = int/float, hi priority = column

  auto words = utils::split_words(format);
  if ((int) words.size() <  size_one)
    error->all(FLERR,"Dump_modify format line is too short");

  int i=0;
  for (const auto &word : words) {
    delete[] vformat[i];

    if (format_column_user[i])
      vformat[i] = utils::strdup(std::string(format_column_user[i]) + " ");
    else if (vtype[i] == Dump::INT && format_int_user)
      vformat[i] = utils::strdup(std::string(format_int_user) + " ");
    else if (vtype[i] == Dump::DOUBLE && format_float_user)
      vformat[i] = utils::strdup(std::string(format_float_user) + " ");
    else if (vtype[i] == Dump::BIGINT && format_bigint_user)
      vformat[i] = utils::strdup(std::string(format_bigint_user) + " ");
    else vformat[i] = utils::strdup(word + " ");
    ++i;
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  if (buffer_flag == 1) write_choice = &DumpLocal::write_string;
  else write_choice = &DumpLocal::write_lines;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  for (i = 0; i < ncompute; i++) {
    int icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump local compute ID");
    compute[i] = modify->compute[icompute];
  }

  for (i = 0; i < nfix; i++) {
    int ifix = modify->find_fix(id_fix[i]);
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
    delete[] label;
    label = utils::strdup(arg[1]);
    return 2;
  } else if (strcmp(arg[0],"format") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");

    if (strcmp(arg[1],"none") == 0) {
      // just clear format_column_user allocated by this dump child class
      for (int i = 0; i < nfield; i++) {
        delete[] format_column_user[i];
        format_column_user[i] = nullptr;
      }
      return 2;
    } else if (strcmp(arg[1],"int") == 0) {
      delete[] format_int_user;
      format_int_user = utils::strdup(arg[2]);
      delete[] format_bigint_user;
      int n = strlen(format_int_user) + 8;
      format_bigint_user = new char[n];
      // replace "d" in format_int_user with bigint format specifier
      // use of &str[1] removes leading '%' from BIGINT_FORMAT string
      char *ptr = strchr(format_int_user,'d');
      if (ptr == nullptr)
        error->all(FLERR,
                   "Dump_modify int format does not contain d character");
      char str[8];
      sprintf(str,"%s",BIGINT_FORMAT);
      *ptr = '\0';
      sprintf(format_bigint_user,"%s%s%s",format_int_user,&str[1],ptr+1);
      *ptr = 'd';

    } else if (strcmp(arg[1],"float") == 0) {
      delete[] format_float_user;
      format_float_user = utils::strdup(arg[2]);

    } else {
      int i = utils::inumeric(FLERR,arg[1],false,lmp) - 1;
      if (i < 0 || i >= nfield)
        error->all(FLERR,"Illegal dump_modify command");
      if (format_column_user[i]) delete[] format_column_user[i];
      format_column_user[i] = utils::strdup(arg[2]);
    }
    return 3;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_header(bigint ndump)
{
  if (me == 0) {
    if (unit_flag && !unit_count) {
      ++unit_count;
      fprintf(fp,"ITEM: UNITS\n%s\n",update->unit_style);
    }
    if (time_flag) fprintf(fp,"ITEM: TIME\n%.16g\n",compute_time());

    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF %s\n",label);
    fprintf(fp,BIGINT_FORMAT "\n",ndump);
    if (domain->triclinic) {
      fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxxlo,boxxhi,boxxy);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxylo,boxyhi,boxxz);
      fprintf(fp,"%-1.16e %-1.16e %-1.16e\n",boxzlo,boxzhi,boxyz);
    } else {
      fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
      fprintf(fp,"%-1.16e %-1.16e\n",boxxlo,boxxhi);
      fprintf(fp,"%-1.16e %-1.16e\n",boxylo,boxyhi);
      fprintf(fp,"%-1.16e %-1.16e\n",boxzlo,boxzhi);
    }
    fprintf(fp,"ITEM: %s %s\n",label,columns);
  }
}

/* ---------------------------------------------------------------------- */

int DumpLocal::count()
{
  int i;

  // invoke Computes for local quantities
  // only if within a run or minimize
  // else require that computes are current
  // this prevents a compute from being invoked by the WriteDump class

  if (ncompute) {
    if (update->whichflag == 0) {
      for (i = 0; i < ncompute; i++)
        if (compute[i]->invoked_local != update->ntimestep)
          error->all(FLERR,"Compute used in dump between runs is not current");
    } else {
      for (i = 0; i < ncompute; i++) {
        if (!(compute[i]->invoked_flag & Compute::INVOKED_LOCAL)) {
          compute[i]->compute_local();
          compute[i]->invoked_flag |= Compute::INVOKED_LOCAL;
        }
      }
    }
  }

  // nmine = # of local values I contribute
  // must be consistent for all input fields

  nmine = -1;

  for (i = 0; i < ncompute; i++) {
    if (nmine < 0) nmine = compute[i]->size_local_rows;
    else if (nmine != compute[i]->size_local_rows)
      error->one(FLERR,
                 "Dump local count is not consistent across input fields");
  }

  for (i = 0; i < nfix; i++) {
    if (nmine < 0) nmine = fix[i]->size_local_rows;
    else if (nmine != fix[i]->size_local_rows)
      error->one(FLERR,
                 "Dump local count is not consistent across input fields");
  }

  return nmine;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::pack(tagint * /*dummy*/)
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpLocal::convert_string(int n, double *mybuf)
{
  int i,j;

  int offset = 0;
  int m = 0;
  for (i = 0; i < n; i++) {
    if (offset + size_one*ONEFIELD > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    for (j = 0; j < size_one; j++) {
      if (vtype[j] == Dump::INT)
        offset += sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE)
        offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
      else if (vtype[j] == Dump::BIGINT)
        offset += sprintf(&sbuf[offset],vformat[j],static_cast<bigint> (mybuf[m]));
      else
        offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
      m++;
    }
    offset += sprintf(&sbuf[offset],"\n");
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_string(int n, double *mybuf)
{
  if (mybuf)
    fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpLocal::write_lines(int n, double *mybuf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == Dump::INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE) fprintf(fp,vformat[j],mybuf[m]);
      else if (vtype[j] == Dump::BIGINT) fprintf(fp,vformat[j],static_cast<bigint>(mybuf[m]));
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

  for (int iarg = 0; iarg < narg; iarg++) {

    if (strcmp(arg[iarg],"index") == 0) {
      pack_choice[iarg] = &DumpLocal::pack_index;
      vtype[iarg] = Dump::INT;

    } else {
      int n;
      ArgInfo argi(arg[iarg],ArgInfo::COMPUTE|ArgInfo::FIX);
      computefixflag = 1;
      vtype[iarg] = Dump::DOUBLE;
      argindex[iarg] = argi.get_index1();

      switch (argi.get_type()) {

        // compute value = c_ID
        // if no trailing [], then arg is set to 0, else arg is int between []

      case ArgInfo::COMPUTE:
        pack_choice[iarg] = &DumpLocal::pack_compute;

        n = modify->find_compute(argi.get_name());
        if (n < 0) error->all(FLERR,"Could not find dump local compute ID");
        if (modify->compute[n]->local_flag == 0)
          error->all(FLERR,"Dump local compute does not compute local info");
        if (argi.get_dim() == 0 && modify->compute[n]->size_local_cols > 0)
          error->all(FLERR,"Dump local compute does not calculate local vector");
        if (argi.get_index1() > 0 && modify->compute[n]->size_local_cols == 0)
          error->all(FLERR,"Dump local compute does not calculate local array");
        if (argi.get_index1() > 0 &&
            argi.get_index1() > modify->compute[n]->size_local_cols)
          error->all(FLERR,"Dump local compute vector is accessed out-of-range");

        field2index[iarg] = add_compute(argi.get_name());
        break;

        // fix value = f_ID
        // if no trailing [], then arg is set to 0, else arg is between []

      case ArgInfo::FIX:
        pack_choice[iarg] = &DumpLocal::pack_fix;

        n = modify->find_fix(argi.get_name());
        if (n < 0) error->all(FLERR,"Could not find dump local fix ID");
        if (modify->fix[n]->local_flag == 0)
          error->all(FLERR,"Dump local fix does not compute local info");
        if (argi.get_dim() == 0 && modify->fix[n]->size_local_cols > 0)
          error->all(FLERR,"Dump local fix does not compute local vector");
        if (argi.get_index1() > 0 && modify->fix[n]->size_local_cols == 0)
          error->all(FLERR,"Dump local fix does not compute local array");
        if (argi.get_index1() > 0 &&
            argi.get_index1() > modify->fix[n]->size_local_cols)
          error->all(FLERR,"Dump local fix vector is accessed out-of-range");

        field2index[iarg] = add_fix(argi.get_name());
        break;

      case ArgInfo::NONE:       // fallthrough
      case ArgInfo::UNKNOWN:    // fallthrough
      default:
         error->all(FLERR,"Invalid attribute in dump local command");
         break;
      }
    }
  }

  if (computefixflag == 0)
    error->all(FLERR,"Dump local attributes contain no compute or fix");
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpLocal::add_compute(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete[] compute;
  compute = new Compute*[ncompute+1];

  id_compute[ncompute] = utils::strdup(id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpLocal::add_fix(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete[] fix;
  fix = new Fix*[nfix+1];

  id_fix[nfix] = utils::strdup(id);
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
