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
#include "dump_cac_xyz.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

#define ONELINE 128
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCACXYZ::DumpCACXYZ(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  typenames(NULL)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xyz command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz filename");

  size_one = 5;

  buffer_allow = 1;
  buffer_flag = 1;
  sort_flag = 0;
  sortcol = 0;

  if (format_default) delete [] format_default;

  char *str = (char *) "%s %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  ntypes = atom->ntypes;
  typenames = NULL;
}

/* ---------------------------------------------------------------------- */

DumpCACXYZ::~DumpCACXYZ()
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

void DumpCACXYZ::init_style()
{
  //check if CAC atom style is defined
  if(!atom->CAC_flag)
  error->all(FLERR, "CAC dump styles require a CAC atom style");
  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

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

  if (buffer_flag == 1) write_choice = &DumpCACXYZ::write_string;
  else write_choice = &DumpCACXYZ::write_lines;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/*------------------------------------------------------------------------*/
int DumpCACXYZ::count()
{
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	int *poly_count = atom->poly_count;


	int **element_scale = atom->element_scale;

	int m = 0;
	for (int i = 0; i < nlocal; i++)
	{

			if (mask[i] & groupbit) m = m + element_scale[i][0]*element_scale[i][1]*element_scale[i][2]*poly_count[i];

	}
	return m;
}
/* ---------------------------------------------------------------------- */

int DumpCACXYZ::modify_param(int narg, char **arg)
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

/* ---------------------------------------------------------------------- */

void DumpCACXYZ::write_header(bigint n)
{
  if (me == 0) {
    fprintf(fp,BIGINT_FORMAT "\n",n);
    fprintf(fp,"Atoms. Timestep: " BIGINT_FORMAT "\n",update->ntimestep);
  }
}
//-------------------------------------------------------------------------

double DumpCACXYZ::shape_function(double s, double t, double w, int flag, int index){
double shape_function=0;
if(flag==2){

    if(index==1){
    shape_function=(1-s)*(1-t)*(1-w)/8;
    }
    else if(index==2){
    shape_function=(1+s)*(1-t)*(1-w)/8;
    }
    else if(index==3){
    shape_function=(1+s)*(1+t)*(1-w)/8;
    }
    else if(index==4){
    shape_function=(1-s)*(1+t)*(1-w)/8;
    }
    else if(index==5){
    shape_function=(1-s)*(1-t)*(1+w)/8;
    }
    else if(index==6){
    shape_function=(1+s)*(1-t)*(1+w)/8;
    }
    else if(index==7){
    shape_function=(1+s)*(1+t)*(1+w)/8;
    }
    else if(index==8){
    shape_function=(1-s)*(1+t)*(1+w)/8;
    }


}
return shape_function;
}

/* ---------------------------------------------------------------------- */

void DumpCACXYZ::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double xmap[3];
  double unit_cell[3];
  double unit_cell_mapped[3];
  double ***current_nodal_positions;
  double shape_func;
  int nodes_per_element;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int nlocal = atom->nlocal;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  double ****nodal_positions = atom->nodal_positions;
  int *periodicity = domain->periodicity;
  double *prd = domain->prd;
  
  //int maptag=1;
  m = n = 0;
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      if(element_type[i]==0){
      buf[m++] = tag[i];
      buf[m++] = node_types[i][0];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      }
      else if(element_type[i]==1){
      unit_cell_mapped[0] = 2 / double(element_scale[i][0]);
			unit_cell_mapped[1] = 2 / double(element_scale[i][1]);
			unit_cell_mapped[2] = 2 / double(element_scale[i][2]);
      current_nodal_positions=nodal_positions[i];
      nodes_per_element=8;

      
        for (int e1 = 0; e1 < element_scale[i][0]; e1++) {
          for (int e2 = 0; e2 < element_scale[i][1]; e2++) {
            for (int e3 = 0; e3 < element_scale[i][2]; e3++) {
              unit_cell[0]=unit_cell_mapped[0]/2+e1*unit_cell_mapped[0]-1;
              unit_cell[1]=unit_cell_mapped[1]/2+e2*unit_cell_mapped[1]-1;
              unit_cell[2]=unit_cell_mapped[2]/2+e3*unit_cell_mapped[2]-1;
              for (int polyscan = 0; polyscan < poly_count[i]; polyscan++) {
                xmap[0]=0;
                xmap[1]=0;
                xmap[2]=0;
                for (int kk = 0; kk < nodes_per_element; kk++) {
            	    shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kk + 1);
				          xmap[0] += current_nodal_positions[polyscan][kk][0] * shape_func;
				          xmap[1] += current_nodal_positions[polyscan][kk][1] * shape_func;
				          xmap[2] += current_nodal_positions[polyscan][kk][2] * shape_func;
                }
              //test if mapped particle is in box and remap otherwise
              if(periodicity[0]){
                if(xmap[0]>boxhi[0]) xmap[0]-=prd[0];
                if(xmap[0]<boxlo[0]) xmap[0]+=prd[0];
              }
              if(periodicity[1]){
                if(xmap[1]>boxhi[1]) xmap[1]-=prd[1];
                if(xmap[1]<boxlo[1]) xmap[1]+=prd[1];
              }
              if(periodicity[2]){
                if(xmap[2]>boxhi[2]) xmap[2]-=prd[2];
                if(xmap[2]<boxlo[2]) xmap[2]+=prd[2];
              }   
              buf[m++] = tag[i];
              buf[m++] = node_types[i][polyscan];
              buf[m++] = xmap[0];
              buf[m++] = xmap[1];
              buf[m++] = xmap[2];
            }
          }
		    }
      }
      }
    if (ids) ids[n++] = tag[i];
    }
  }
}


/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCACXYZ::convert_string(int n, double *mybuf)
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
                      typenames[static_cast<int> (mybuf[m+1])],
                      mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpCACXYZ::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCACXYZ::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCACXYZ::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            typenames[static_cast<int> (mybuf[m+1])],
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}
