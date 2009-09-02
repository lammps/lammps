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

/* ----------------------------------------------------------------------
   Contributing author: Liang Wan (Chinese Academy of Sciences)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "dump_cfg.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "fix.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{INT,DOUBLE};  // same as in dump_custom.cpp

/* ---------------------------------------------------------------------- */

DumpCFG::DumpCFG(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  if (narg < 10 ||
      strcmp(arg[5],"id") != 0 || strcmp(arg[6],"type") != 0 ||
      strcmp(arg[7],"xs") != 0 || strcmp(arg[8],"ys") != 0 ||
      strcmp(arg[9],"zs") != 0)
    error->all("Dump cfg arguments must start with 'id type xs ys zs'");

  ntypes = atom->ntypes;
  typenames = NULL;

  // arrays for data rearrangement

  recvcounts = new int[nprocs];
  displs = new int[nprocs];
  tags = NULL;
  rbuf = NULL;
  nchosen = nlines = 0;

  // setup auxiliary property name strings
  // convert 'X_ID[m]' (X=c,f,v) to 'ID_m'

  if (narg > 10) auxname = new char*[narg-10];
  else auxname = NULL;

  int i = 0;
  for (int iarg = 10; iarg < narg; iarg++, i++) {
    if (strncmp(arg[iarg],"c_",2) == 0 ||
	strncmp(arg[iarg],"f_",2) == 0 ||
	strncmp(arg[iarg],"v_",2) == 0) {
      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump cfg command");
	*ptr = '\0';
	*(ptr+2) = '\0';
	auxname[i] = new char[strlen(suffix) + 3];
	strcpy(auxname[i],suffix);
	strcat(auxname[i],"_");
	strcat(auxname[i],ptr+1);
      } else {
	auxname[i] = new char[strlen(suffix) + 1];
	strcpy(auxname[i],suffix);
      }

      delete [] suffix;

    } else {
      auxname[i] = new char[strlen(arg[iarg]) + 1];
      strcpy(auxname[i],arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

DumpCFG::~DumpCFG()
{
  if (typenames) {
    for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
    delete [] typenames;
  }

  delete [] recvcounts;
  delete [] displs;

  if (tags) memory->sfree(tags);
  if (rbuf) memory->destroy_2d_double_array(rbuf);

  if (auxname) {
    for (int i = 0; i < nfield-5; i++) delete [] auxname[i];
    delete [] auxname;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCFG::init()
{
  if (multifile == 0)
    error->all("Dump in CFG format requires one snapshot per file");

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[3];
      strcpy(typenames[itype],"C");
    }
    if (comm->me == 0)
      error->warning("All element names have been set to 'C' for dump cfg");
  }

  // setup format strings

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
    if (ptr == NULL) error->all("Dump cfg user format string error");
    delete [] vformat[i];
    vformat[i] = new char[strlen(ptr) + 2];
    strcpy(vformat[i],ptr);
    vformat[i] = strcat(vformat[i]," ");
  }
  if (strtok(NULL," \0"))
    error->all("Dump cfg user format string error");

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all("Could not find dump cfg compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all("Could not find dump cfg fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all("Dump cfg and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) error->all("Could not find dump cfg variable name");
    variable[i] = ivariable;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCFG::write_header(int n)
{
  if (me == 0 || multiproc) {
    fprintf(fp,"Number of particles = %d\n", n);
    fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
    fprintf(fp,"H0(1,1) = %g A\n",domain->xprd);
    fprintf(fp,"H0(1,2) = 0 A \n");
    fprintf(fp,"H0(1,3) = 0 A \n");
    fprintf(fp,"H0(2,1) = %g A \n",domain->xy);
    fprintf(fp,"H0(2,2) = %g A\n",domain->yprd);
    fprintf(fp,"H0(2,3) = 0 A \n");
    fprintf(fp,"H0(3,1) = %g A \n",domain->xz);
    fprintf(fp,"H0(3,2) = %g A \n",domain->yz);
    fprintf(fp,"H0(3,3) = %g A\n",domain->zprd);
    fprintf(fp,".NO_VELOCITY.\n");
    fprintf(fp,"entry_count = %d\n",nfield-2);
    for (int i = 0; i < nfield-5; i++)
      fprintf(fp,"auxiliary[%d] = %s\n",i,auxname[i]);
  }

  // calculate total # of data lines to be written on a writing proc

  if (multiproc) nchosen = nmine;
  else MPI_Reduce(&nmine,&nchosen,1,MPI_INT,MPI_SUM,0,world);

  // allocate memory needed for data rearrangement on writing proc(s)

  if (multiproc || me == 0) {
    if (rbuf) memory->destroy_2d_double_array(rbuf);
    rbuf = memory->create_2d_double_array(nchosen,size_one,"dump:rbuf");
  }

  // create a sorted list of atom IDs on writing proc(s) if necessary

  if (sort_flag) create_sorted_tags();
}

/* ----------------------------------------------------------------------
   write data lines to file in a block-by-block style
   write head of block (mass & element name) only if has atoms of the type
------------------------------------------------------------------------- */

void DumpCFG::write_data(int n, double *buf)
{
  int i,j,m,itype;
  int tag_i,index;
  double *mass = atom->mass;

  // if sort flag is off, transfer data in buf to rbuf directly
  // if write by proc 0, transfer chunk by chunk

  if (sort_flag == 0) {
    for (i = 0, m = 0; i < n; i++) {
      for (j = 0; j < size_one; j++)
	rbuf[nlines][j] = buf[m++];
      nlines++;
    }

    //  write data lines in rbuf to file after transfer is done

    if (nlines == nchosen) {
      for (itype = 1; itype <= ntypes; itype++) {
	for (i = 0; i < nchosen; i++)
	  if (rbuf[i][1] == itype) break;
	if (i < nchosen) {
	  fprintf(fp,"%g\n",mass[itype]);
	  fprintf(fp,"%s\n",typenames[itype]);
	  for (; i < nchosen; i++) {
	    if (rbuf[i][1] == itype) {
	      for (j = 2; j < size_one; j++) {
		if (vtype[j] == INT)
		  fprintf(fp,vformat[j],static_cast<int> (rbuf[i][j]));
		else fprintf(fp,vformat[j],rbuf[i][j]);
	      }
	      fprintf(fp,"\n");
	    }
	  }
	}
      }
      nlines = 0;
    }

  // if sort flag is on, transfer data lines in buf to rbuf &
  //   rearrange them in a sorted order
  // if write by proc 0, transfer chunk by chunk

  } else {

    // sort data lines:
    // find the index of the atom ID of each data line within the
    //   sorted atom IDs array
    // transfer the data line to rbuf with its location according
    //   to this index

    for (i = 0, m = 0; i < n; i++) {
      tag_i = static_cast<int>(buf[m]);
      for (j = 0; j < nchosen; j++) {
	if (tag_i == tags[j]) {
	  index = j;
	  break;
	}
      }
      for (j = 0; j < size_one; j++)
	rbuf[index][j] = buf[m++];
    }

    //  write data lines in rbuf to file after transfer is done

    nlines += n;
    if (nlines == nchosen) {
      for (itype = 1; itype <= ntypes; itype++) {
	for (i = 0; i < nchosen; i++)
	  if (rbuf[i][1] == itype) break;
	if (i < nchosen) {
	  fprintf(fp,"%g\n",mass[itype]);
	  fprintf(fp,"%s\n",typenames[itype]);
	  for (; i < nchosen; i++) {
	    if (rbuf[i][1] == itype) {
	      for (j = 2; j < size_one; j++) {
		if (vtype[j] == INT)
		  fprintf(fp,vformat[j],static_cast<int>(rbuf[i][j]));
		else fprintf(fp,vformat[j],rbuf[i][j]);
	      }
	      fprintf(fp,"\n");
	    }
	  }
	}
      }
      nlines = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   create a sorted list of atom IDs on writing proc(s)
------------------------------------------------------------------------- */

void DumpCFG::create_sorted_tags()
{
  int index;
  int *mytags;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  index = 0;

  // if write by multiproc, each proc has its own atom IDs collected

  if (multiproc) {
    if (tags) memory->sfree(tags);
    tags = (int *) memory->smalloc(nchosen*sizeof(int),"dump:tags");
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) tags[index++] = tag[i];

  // if write by proc 0, gather atom IDs of all the chosen atoms on proc 0

  } else {
    mytags = (int *) memory->smalloc(nmine*sizeof(int),"dump:mytags");
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) mytags[index++] = tag[i];
    MPI_Gather(&nmine,1,MPI_INT,recvcounts,1,MPI_INT,0,world);
    if (me == 0) {
      if (tags) memory->sfree(tags);
      tags = (int *) memory->smalloc(nchosen*sizeof(int),"dump:tags");
      displs[0] = 0;
      for (int iproc = 1; iproc < nprocs; iproc++)
	displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
    }
    MPI_Gatherv(mytags,nmine,MPI_INT,tags,recvcounts,displs,MPI_INT,0,world);
    memory->sfree(mytags);
  }

  if (multiproc || me == 0) qsort(tags,nchosen,sizeof(int),tag_compare);
}

/* ----------------------------------------------------------------------
   compare function for qsort
------------------------------------------------------------------------- */

int DumpCFG::tag_compare(const void *itag, const void *jtag)
{
  int tag_i = *((int *) itag);
  int tag_j = *((int *) jtag);

  if (tag_i < tag_j) return -1;
  else if (tag_i > tag_j) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

int DumpCFG::modify_param2(int narg, char **arg)
{
  if (strcmp(arg[0],"element") == 0) {
    if (narg != ntypes+1)
      error->all("Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[3];
      if (strlen(arg[itype]) >= 3)
	error->all("Illegal chemical element names");
      strcpy(typenames[itype],arg[itype]);
    }
    return ntypes+1;

  } else return 0;
}
