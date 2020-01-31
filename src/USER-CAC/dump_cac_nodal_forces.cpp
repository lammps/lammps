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

#include <cstring>
#include "dump_cac_nodal_forces.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

#define ONELINE 128
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCACNodalForces::DumpCACNodalForces(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xyz command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz filename");

  size_one = 6;

  buffer_allow = 1;
  buffer_flag = 1;
  sort_flag = 0;
  sortcol = 0;

  if (format_default) delete [] format_default;

  char *str = (char *) "%d %d %d %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  ntypes = atom->ntypes;
  typenames = NULL;
}

/* ---------------------------------------------------------------------- */

DumpCACNodalForces::~DumpCACNodalForces()
{
  delete[] format_default;
  format_default = NULL;

  if (typenames) {
    for (int i = 1; i <= ntypes; i++)
      delete [] typenames[i];
    delete [] typenames;
    typenames = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::init_style()
{
  //check if CAC atom style is defined
  if(!atom->CAC_flag)
  error->all(FLERR, "CAC dump styles require a CAC atom style");

  //check if sorting was used
  if(sort_flag)
  error->all(FLERR, "CAC dump styles cannot currently be sorted");

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");
  nodes_per_element = atom->nodes_per_element;
  maxpoly = atom->maxpoly;
  // initialize typenames array to be backward compatible by default
  // a 32-bit int can be maximally 10 digits plus sign

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[12];
      sprintf(typenames[itype],"%d",itype);
    }
  }

  // setup function ptr

  if (buffer_flag == 1) write_choice = &DumpCACNodalForces::write_string;
  else write_choice = &DumpCACNodalForces::write_lines;

  // open single file, one time only

  if (multifile == 0) openfile();
  ptimestep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int DumpCACNodalForces::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR, "Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++)
        delete [] typenames[i];

      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 1;
      typenames[itype] = new char[n];
      strcpy(typenames[itype],arg[itype]);
    }

    return ntypes+1;
  }

  return 0;
}


/*------------------------------------------------------------------------*/
int DumpCACNodalForces::count()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *element_type= atom->element_type;
  int *poly_count = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int m = 0;

  //compute number of nodes in total system
  int local_node_count=0;
  total_node_count=0;
  int local_element_count=0;
  total_element_count=0;

  for (int i=0; i<atom->nlocal; i++){
    if (mask[i] & groupbit){
    local_node_count+=nodes_per_element_list[element_type[i]];
    local_element_count++;
    }
  }
  MPI_Allreduce(&local_node_count,&total_node_count,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&local_element_count,&total_element_count,1,MPI_INT,MPI_SUM,world);


  for (int i = 0; i < nlocal; i++)
  {
    if (update->ntimestep - ptimestep == 0) {
      if (mask[i] & groupbit) m = m + nodes_per_element_list[element_type[i]]*poly_count[i] + 1;
    }
    else {
      if (mask[i] & groupbit) m = m + nodes_per_element_list[element_type[i]]*poly_count[i] + 1;
    }
  }
  return m;
}
/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::write_header(bigint n)
{

  if (me == 0) {
  fprintf(fp, " t= " BIGINT_FORMAT " n= " BIGINT_FORMAT
  " e= " BIGINT_FORMAT " Q4 " "\n",
  update->ntimestep, (bigint)total_node_count, total_element_count);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double ****nodal_forces = atom->nodal_forces;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int nlocal = atom->nlocal;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] =double( tag[i]);
      buf[m++] = double( element_type[i]);
      buf[m++] = double (poly_count[i]);
      buf[m++] = double(element_scale[i][0]);
      buf[m++] = double(element_scale[i][1]);
      buf[m++] = double(element_scale[i][2]);

    for (int k = 0; k < poly_count[i]; k++) {
      for (int j = 0; j < nodes_per_element_list[element_type[i]]; j++) {
        buf[m++] = double(j + 1);
        buf[m++] = double(k + 1);
        buf[m++] = double(node_types[i][k]);
        buf[m++] = nodal_forces[i][k][j][0];
        buf[m++] = nodal_forces[i][k][j][1];
        buf[m++] = nodal_forces[i][k][j][2];
      }
    }
    }
  }
}


/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCACNodalForces::convert_string(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
    static_cast<tagint> (mybuf[m]),
    static_cast<tagint>(mybuf[m+1]), static_cast<tagint>(mybuf[m+2]),
    mybuf[m+3],mybuf[m+4],mybuf[m+5]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalForces::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
      typenames[static_cast<int> (mybuf[m+1])],
      mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }


}
