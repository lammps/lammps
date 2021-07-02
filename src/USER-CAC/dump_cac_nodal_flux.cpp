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

#include "dump_cac_nodal_flux.h"
#include <cstring>
#include <cmath>
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "force.h"

using namespace LAMMPS_NS;

#define ONELINE 2048
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCACNodalFlux::DumpCACNodalFlux(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg < 9 || narg > 12) error->all(FLERR,"Illegal dump cac/fluxes command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump cac/fluxes filename");

  if (strcmp(arg[5], "nodal") == 0) cac_flux_flag = 1;
  else if (strcmp(arg[5], "quadrature") == 0) cac_flux_flag = 2;
  else error->all(FLERR, "Unexpected argument in dump cac/flux invocation");
  
  //set atom class flag that can be used by other classes
  atom->cac_flux_flag = cac_flux_flag;

  box_size[0]=0;
  box_size[1]=0;
  box_size[2]=0;
  box_center[0]=0;
  box_center[1]=0;
  box_center[2]=0;

  box_size[0]=utils::numeric(FLERR,arg[6],false,lmp);
  box_size[1]=utils::numeric(FLERR,arg[7],false,lmp);
  box_size[2]=utils::numeric(FLERR,arg[8],false,lmp);
  
  if(narg>9){
    if(narg!=12)
    error->all(FLERR,"Illegal dump cac/fluxes command; if specifying box centers specify all three dims.");
    box_center[0]=utils::numeric(FLERR,arg[9],false,lmp);
    box_center[1]=utils::numeric(FLERR,arg[10],false,lmp);
    box_center[2]=utils::numeric(FLERR,arg[11],false,lmp);
  }

  //test box sizes for feasibility
  if(box_size[0]<=0||box_size[1]<=0||box_size[2]<=0)
  error->all(FLERR,"box sizes for fix surface/fluxes/atom must be greater than zero");

  //make sure box_center displacement is smaller than half the box size (so that the atom stays inside the box)
  if((box_center[0]>0&&box_size[0]<=box_center[0])||(box_center[0]<0&&-box_size[0]>=box_center[0]))
  error->all(FLERR,"box center cannot be larger than or equal to half the box size (atom must remain inside the box)");
  if((box_center[1]>0&&box_size[1]<=box_center[1])||(box_center[1]<0&&-box_size[1]>=box_center[1]))
  error->all(FLERR,"box center cannot be larger than or equal to half the box size (atom must remain inside the box)");
  if((box_center[2]>0&&box_size[2]<=box_center[2])||(box_center[2]<0&&-box_size[2]>=box_center[2]))
  error->all(FLERR,"box center cannot be larger than or equal to half the box size (atom must remain inside the box)");

  size_one = 27;

  buffer_allow = 1;
  buffer_flag = 1;
  sort_flag = 0;
  sortcol = 0;

  if (format_default) delete [] format_default;

  char *str = (char *) "%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  ntypes = atom->ntypes;
  typenames = NULL;
}

/* ---------------------------------------------------------------------- */

DumpCACNodalFlux::~DumpCACNodalFlux()
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

void DumpCACNodalFlux::init_style()
{
  //check if CAC atom style is defined
  if(!atom->CAC_flag)
  error->all(FLERR, "CAC dump styles require a CAC atom style");
  
  //set box parameters that can be used by the cac pair style
  atom->box_size[0] = box_size[0];
  atom->box_size[1] = box_size[1];
  atom->box_size[2] = box_size[2];
  atom->box_center[0] = box_center[0];
  atom->box_center[1] = box_center[1];
  atom->box_center[2] = box_center[2];

  //compute additional cutoff that will be needed to ensure quadrature points have
  //  the necessary neighbor list
  double center_distance = sqrt(box_center[0]*box_center[0]+box_center[1]*box_center[1]
  +box_center[2]*box_center[2]);
  double box_distance = 0.5*sqrt(box_size[0]*box_size[0]+box_size[1]*box_size[1]
  +box_size[2]*box_size[2]);
  atom->cut_add = center_distance + box_distance;

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

  if (buffer_flag == 1) write_choice = &DumpCACNodalFlux::write_string;
  else write_choice = &DumpCACNodalFlux::write_lines;

  // open single file, one time only

  if (multifile == 0) openfile();
  ptimestep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int DumpCACNodalFlux::modify_param(int narg, char **arg)
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
int DumpCACNodalFlux::count()
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

void DumpCACNodalFlux::write_header(bigint n)
{

  if (me == 0) {
  fprintf(fp, " t= " BIGINT_FORMAT " n= " BIGINT_FORMAT
  " e= " BIGINT_FORMAT " Q4 " "\n",
  update->ntimestep, (bigint)total_node_count, total_element_count);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalFlux::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double ****nodal_fluxes = atom->nodal_fluxes;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int nlocal = atom->nlocal;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  double nktv2p = force->nktv2p;
  
  m = n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      buf[m++] =double( tag[i]);
      buf[m++] = double( element_type[i]);
      buf[m++] = double (poly_count[i]);
      buf[m++] = double(element_scale[i][0]);
      buf[m++] = double(element_scale[i][1]);
      buf[m++] = double(element_scale[i][2]);
      //pad remainder of this line with 0s due to current dump constraints
      for(int ipad=0; ipad < size_one-6; ipad++) buf[m++] = 0;

    for (int k = 0; k < poly_count[i]; k++) {
      for (int j = 0; j < nodes_per_element_list[element_type[i]]; j++) {
        buf[m++] = double(j + 1);
        buf[m++] = double(k + 1);
        buf[m++] = double(node_types[i][k]);
        for(int ibuf=0; ibuf < size_one-3; ibuf++)
        if(ibuf%4!=0)
          buf[m++] = nktv2p*nodal_fluxes[i][k][j][ibuf];
        else
          buf[m++] = nodal_fluxes[i][k][j][ibuf];
      }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCACNodalFlux::convert_string(int n, double *mybuf)
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
      mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9],mybuf[m+10],mybuf[m+11],
      mybuf[m+12],mybuf[m+13],mybuf[m+14],mybuf[m+15],mybuf[m+16],mybuf[m+17],mybuf[m+18],
      mybuf[m+19],mybuf[m+20],mybuf[m+21],mybuf[m+22],mybuf[m+23],mybuf[m+24],mybuf[m+25],mybuf[m+26]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalFlux::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalFlux::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCACNodalFlux::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
      typenames[static_cast<int> (mybuf[m+1])],
      mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9],mybuf[m+10],mybuf[m+11],
      mybuf[m+12],mybuf[m+13],mybuf[m+14],mybuf[m+15],mybuf[m+16],mybuf[m+17],mybuf[m+18],
      mybuf[m+19],mybuf[m+20],mybuf[m+21],mybuf[m+22],mybuf[m+23],mybuf[m+24],mybuf[m+25],mybuf[m+26]);
    m += size_one;
  }


}
