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
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#include "inthash.h"

using namespace LAMMPS_NS;

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

  size_one = 0;
  maxbuf = 0;
  buf = NULL;

  idmap = NULL;
  size_one_id = sizeof(inthash_node_t);

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

  // delete hashtable, if it exists
  release_idmap();

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

/* ----------------------------------------------------------------------
   build hash table for sorted output
------------------------------------------------------------------------- */
void Dump::build_idmap(int natoms)
{
  /* nme:    number of atoms in group on this MPI task
   * nmax:   max number of atoms in group across all MPI tasks
   * nlocal: all local atoms */
  int nmax,nme,nlocal;
  int *mask  = atom->mask;
  int *tag  = atom->tag;
  nlocal = atom->nlocal;

  int i,j;

  nme=0;
  for (i=0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nme;

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  int *comm_buf = new int[nmax];
  int comm_size = nmax*sizeof(int);

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;

  if (me == 0) {
    if (screen)
      fprintf(screen,"Preparing index map for sorted dump of %d atoms\n",natoms);
    if (logfile)
      fprintf(logfile,"Preparing index map for sorted dump of %d atoms\n",natoms);

    /* initialize and build hashtable. */
    inthash_t *hashtable=new inthash_t;
    inthash_init(hashtable, natoms);
    idmap = (void *)hashtable;
    int *taglist = new int[natoms];
    int numtag=0; /* counter to map atoms to a 0-based consecutive index */
    
    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        taglist[numtag] = tag[i];
        ++numtag;
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      MPI_Irecv(comm_buf, comm_size, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= sizeof(int);

      for (j=0; j < ndata; ++j) {
        taglist[numtag] = comm_buf[j];
        ++numtag;
      }
    }

    /* sort list of tags by value to have consistently the
     * same list when running in parallel and build hash table. */
    id_sort(taglist, 0, natoms-1);
    for (i=0; i < natoms; ++i) {
      inthash_insert(hashtable, taglist[i], i);
    }
    delete[] taglist;

  } else {
    nme=0;
    for (i=0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        comm_buf[nme] = tag[i];
        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, &status);
    MPI_Rsend(comm_buf, nme*sizeof(int), MPI_BYTE, 0, 0, world);
  }
  // cleanup.
  delete[] comm_buf;
}

/* ----------------------------------------------------------------------
   delete hashtable 
------------------------------------------------------------------------- */

void Dump::release_idmap(void)
{
  if (idmap) {
    inthash_t *hashtable = (inthash_t *)idmap;
    inthash_destroy(hashtable);
    delete hashtable;
    hashtable = NULL;
  }
}

/* ----------------------------------------------------------------------
   lookup id in hashtable.
------------------------------------------------------------------------- */

int Dump::lookup_id(int tag)
{
  if (idmap) {
    inthash_t *hashtable = (inthash_t *)idmap;
    return inthash_lookup(hashtable, tag);
  } else {
    return tag;
  }
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
  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int nme = count();

  int ntotal,nmax;
  if (multiproc) nmax = nme;
  else {
    MPI_Allreduce(&nme,&ntotal,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  }

  // write timestep header

  if (multiproc) write_header(nme);
  else write_header(ntotal);

  // grow communication buffer if necessary

  if (nmax*size_one > maxbuf) {
    maxbuf = nmax*size_one;
    memory->sfree(buf);
    buf = (double *) memory->smalloc(maxbuf*sizeof(double),"dump:buf");
  }

  // pack my data into buf
  // me_size = # of quantities in buf

  int me_size = pack();

  // multiproc = 1 = each proc writes own data to own file 
  // multiproc = 0 = all procs write to one file thru proc 0
  //   proc 0 pings each proc, receives it's data, writes to file
  //   all other procs wait for ping, send their data to proc 0

  if (multiproc) write_data(me_size/size_one,buf);
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
	} else nlines = me_size/size_one;

	write_data(nlines,buf);
      }
      if (flush_flag) fflush(fp);
      
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,me_size,MPI_DOUBLE,0,0,world);
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
      if (strcmp(arg[iarg+1],"yes") == 0) sort_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) sort_flag = 0;
      else error->all("Illegal dump_modify command");
      iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all("Illegal dump_modify command");
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf
------------------------------------------------------------------------- */

double Dump::memory_usage()
{
  double bytes = maxbuf * sizeof(double);
  return bytes;
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
