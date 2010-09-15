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
#include "stdio.h"
#include "dump.h"
#include "atom.h"
#include "irregular.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// allocate space for static class variable

Dump *Dump::dumpptr;

#define BIG 1.0e20
#define IBIG 2147483647
#define EPSILON 1.0e-6

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Dump::Dump(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  first_flag = 0;
  flush_flag = 1;
  format = NULL;
  format_user = NULL;
  clearstep = 0;
  sort_flag = 0;
  append_flag = 0;

  maxbuf = maxsort = maxproc = 0;
  buf = bufsort = NULL;
  idsort = index = proclist = NULL;
  irregular = NULL;

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  fp = NULL;
  singlefile_opened = 0;
  compressed = 0;
  binary = 0;
  multifile = 0;
  multiproc = 0;

  char *ptr;
  if (ptr = strchr(filename,'%')) {
    multiproc = 1;
    char *extend = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(extend,"%s%d%s",filename,me,ptr+1);
    delete [] filename;
    n = strlen(extend) + 1;
    filename = new char[n];
    strcpy(filename,extend);
    delete [] extend;
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;

  delete [] format;
  delete [] format_default;
  delete [] format_user;

  memory->sfree(buf);
  memory->sfree(bufsort);
  memory->sfree(idsort);
  memory->sfree(index);
  memory->sfree(proclist);
  delete irregular;

  // XTC style sets fp to NULL since it closes file in its destructor

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();

  if (!sort_flag) {
    memory->sfree(buf);
    memory->sfree(bufsort);
    memory->sfree(idsort);
    memory->sfree(index);
    memory->sfree(proclist);
    delete irregular;

    maxbuf = maxsort = maxproc = 0;
    buf = bufsort = NULL;
    idsort = index = proclist = NULL;
    irregular = NULL;
  }

  if (sort_flag && sortcol == 0 && atom->tag_enable == 0)
      error->all("Cannot use dump sort on atom IDs with no atom IDs defined");

  if (sort_flag && sortcol > size_one)
    error->all("Dump sort column is invalid");

  if (sort_flag && nprocs > 1 && irregular == NULL)
    irregular = new Irregular(lmp);
}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc will contribute to dump

  nme = count();

  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int ntotal,nmax;
  if (multiproc) nmax = nme;
  else {
    MPI_Allreduce(&nme,&ntotal,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  }

  // write timestep header

  if (multiproc) write_header(nme);
  else write_header(ntotal);

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nmax > maxbuf) {
    maxbuf = nmax;
    memory->sfree(buf);
    buf = (double *) 
      memory->smalloc(maxbuf*size_one*sizeof(double),"dump:buf");
    if (sort_flag && sortcol == 0) {
      memory->sfree(ids);
      ids = (int *) memory->smalloc(maxbuf*sizeof(int),"dump:ids");
    }
  }

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(NULL);
  if (sort_flag) sort();

  // multiproc = 1 = each proc writes own data to own file 
  // multiproc = 0 = all procs write to one file thru proc 0
  //   proc 0 pings each proc, receives it's data, writes to file
  //   all other procs wait for ping, send their data to proc 0

  if (multiproc) write_data(nme,buf);
  else {
    int tmp,nlines;
    MPI_Status status;
    MPI_Request request;

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
	if (iproc) {
	  MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	  nlines /= size_one;
	} else nlines = nme;

	write_data(nlines,buf);
      }
      if (flush_flag) fflush(fp);
      
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
    }
  }

  // if file per timestep, close file

  if (multifile) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent;
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    sprintf(filecurrent,"%s%d%s",filename,update->ntimestep,ptr+1);
    *ptr = '*';
  }

  // open one file on proc 0 or file on every proc

  if (me == 0 || multiproc) {
    if (compressed) {
#ifdef LAMMPS_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
      fp = popen(gzip,"w");
#else
      error->one("Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    } else if (append_flag) {
      fp = fopen(filecurrent,"a");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == NULL) error->one("Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ----------------------------------------------------------------------
   parallel sort of buf across all procs
   changes nme, reorders datums in buf, grows buf if necessary
------------------------------------------------------------------------- */

void Dump::sort()
{
  int i,iproc;
  double value;

  // if single proc, swap ptrs to buf,ids <-> bufsort,idsort

  if (nprocs == 1) {
    if (nme > maxsort) {
      maxsort = nme;
      memory->sfree(bufsort);
      bufsort = (double *)
	memory->smalloc(maxsort*size_one*sizeof(double),"dump:bufsort");
      memory->sfree(index);
      index = (int *) memory->smalloc(maxsort*sizeof(int),"dump:index");
      if (sortcol == 0) {
	memory->sfree(idsort);
	idsort = (int *) memory->smalloc(maxsort*sizeof(int),"dump:idsort");
      }
    }

    double *dptr = buf;
    buf = bufsort;
    bufsort = dptr;

    if (sortcol == 0) {
      int *iptr = ids;
      ids = idsort;
      idsort = iptr;
    }

  // if multiple procs, exchange datums between procs via irregular
    
  } else {

    // grow proclist if necessary

    if (nme > maxproc) {
      maxproc = nme;
      memory->sfree(proclist);
      proclist = (int *) memory->smalloc(maxproc*sizeof(int),"dump:proclist");
    }
    
    // proclist[i] = which proc Ith datum will be sent to

    if (sortcol == 0) {
      int min = IBIG;
      int max = 0;
      for (i = 0; i < nme; i++) {
	min = MIN(min,ids[i]);
	max = MAX(max,ids[i]);
      }
      int minall,maxall;
      MPI_Allreduce(&min,&minall,1,MPI_INT,MPI_MIN,world);
      MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
      double range = maxall - minall + 0.1;
      for (i = 0; i < nme; i++) {
	iproc = static_cast<int> ((ids[i]-minall)/range * nprocs);
	proclist[i] = iproc;
      }

    } else {
      double min = BIG;
      double max = -BIG;
      for (i = 0; i < nme; i++) {
	value = buf[i*size_one + sortcolm1];
	min = MIN(min,value);
	max = MAX(max,value);
      }
      double minall,maxall;
      MPI_Allreduce(&min,&minall,1,MPI_DOUBLE,MPI_MIN,world);
      MPI_Allreduce(&max,&maxall,1,MPI_DOUBLE,MPI_MAX,world);
      double range = maxall-minall + EPSILON*(maxall-minall);
      if (range == 0.0) range = EPSILON;
      for (i = 0; i < nme; i++) {
	value = buf[i*size_one + sortcolm1];
	iproc = static_cast<int> ((value-minall)/range * nprocs);
	proclist[i] = iproc;
      }
    }

    // create comm plan, grow recv bufs if necessary,
    // exchange datums, destroy plan
    // if sorting on atom IDs, exchange IDs also

    nme = irregular->create_data(nme,proclist);

    if (nme > maxsort) {
      maxsort = nme;
      memory->sfree(bufsort);
      bufsort = (double *) 
	memory->smalloc(maxsort*size_one*sizeof(double),"dump:bufsort");
      memory->sfree(index);
      index = (int *) memory->smalloc(maxsort*sizeof(int),"dump:index");
      if (sortcol == 0) {
	memory->sfree(idsort);
	idsort = (int *) memory->smalloc(maxsort*sizeof(int),"dump:idsort");
      }
    }
    
    irregular->exchange_data((char *) buf,size_one*sizeof(double),
			     (char *) bufsort);
    if (sortcol == 0) 
      irregular->exchange_data((char *) ids,sizeof(int),(char *) idsort);
    irregular->destroy_data();
  }

  // quicksort indices using IDs or buf column as comparator

  dumpptr = this;
  for (i = 0; i < nme; i++) index[i] = i;
  if (sortcol == 0) qsort(index,nme,sizeof(int),idcompare);
  else qsort(index,nme,sizeof(int),bufcompare);

  // copy data from bufsort to buf using index

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->sfree(buf);
    buf = (double *) 
      memory->smalloc(maxbuf*size_one*sizeof(double),"dump:buf");
  }
    
  int nbytes = size_one*sizeof(double);
  for (i = 0; i < nme; i++)
    memcpy(&buf[i*size_one],&bufsort[index[i]*size_one],nbytes);
}

/* ----------------------------------------------------------------------
   compare two atom IDs in sort data structure
   called via qsort_r in sort() method
   is a static method so access sort data structure via ptr
------------------------------------------------------------------------- */

int Dump::idcompare(const void *pi, const void *pj)
{
  int *idsort = dumpptr->idsort;

  int i = *((int *) pi);
  int j = *((int *) pj);

  if (idsort[i] < idsort[j]) return -1;
  if (idsort[i] > idsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compare two buffer quantities in sort data structure with size_one stride
   called via qsort_r in sort() method
   is a static method so access sort data structure via ptr
------------------------------------------------------------------------- */

int Dump::bufcompare(const void *pi, const void *pj)
{
  double *bufsort = dumpptr->bufsort;
  int size_one = dumpptr->size_one;
  int sortcolm1 = dumpptr->sortcolm1;

  int i = *((int *) pi)*size_one + sortcolm1;
  int j = *((int *) pj)*size_one + sortcolm1;

  if (bufsort[i] < bufsort[j]) return -1;
  if (bufsort[i] > bufsort[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process params common to all dumps here
   if unknown param, call modify_param specific to the dump
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) append_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) append_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      int idump;
      for (idump = 0; idump < output->ndump; idump++)
	if (strcmp(id,output->dump[idump]->id) == 0) break;
      int n;
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	delete [] output->var_dump[idump];
	n = strlen(&arg[iarg+1][2]) + 1;
	output->var_dump[idump] = new char[n];
	strcpy(output->var_dump[idump],&arg[iarg+1][2]);
	n = 0;
      } else {
	n = atoi(arg[iarg+1]);
	if (n <= 0) error->all("Illegal dump_modify command");
      }
      output->every_dump[idump] = n;
      iarg += 2;
    } else if (strcmp(arg[iarg],"first") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) first_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) first_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      delete [] format_user;
      format_user = NULL;
      if (strcmp(arg[iarg+1],"none")) {
	int n = strlen(arg[iarg+1]) + 1;
	format_user = new char[n];
	strcpy(format_user,arg[iarg+1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"sort") == 0) {
      if (iarg+2 > narg) error->all("Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"off") == 0) sort_flag = 0;
      else {
	sort_flag = 1;
	sortcol = atoi(arg[iarg+1]);
	if (sortcol < 0) error->all("Illegal dump_modify command");
	sortcolm1 = sortcol - 1;
      }
      iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all("Illegal dump_modify command");
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double Dump::memory_usage()
{
  double bytes = maxbuf*size_one * sizeof(double);
  if (sort_flag) {
    if (sortcol == 0) bytes += maxbuf * sizeof(int);    // ids
    bytes += maxsort*size_one * sizeof(double);         // bufsort
    bytes += maxsort * sizeof(int);                     // index
    if (sortcol == 0) bytes += maxsort * sizeof(int);   // idsort
    bytes += maxproc * sizeof(int);                     // proclist
    if (irregular) bytes += irregular->memory_usage();
  }
  return bytes;
}
