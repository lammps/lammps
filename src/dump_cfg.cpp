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
   	Memory efficiency improved by Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
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

enum{INT,DOUBLE,STRING,BIGINT};   // same as in DumpCustom

#define UNWRAPEXPAND 10.0
#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpCFG::DumpCFG(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  if (narg < 10 ||
      strcmp(arg[5],"mass") != 0 || strcmp(arg[6],"type") != 0 ||
      (strcmp(arg[7],"xs") != 0 && strcmp(arg[7],"xsu") != 0) ||
      (strcmp(arg[8],"ys") != 0 && strcmp(arg[8],"ysu") != 0) ||
      (strcmp(arg[9],"zs") != 0 && strcmp(arg[9],"zsu") != 0))
    error->all(FLERR,"Dump cfg arguments must start with "
               "'mass type xs ys zs' or 'mass type xsu ysu zsu'");

  if (strcmp(arg[7],"xs") == 0)
    if (strcmp(arg[8],"ysu") == 0 || strcmp(arg[9],"zsu") == 0)
      error->all(FLERR,
                 "Dump cfg arguments can not mix xs|ys|zs with xsu|ysu|zsu");
    else unwrapflag = 0;
  else if (strcmp(arg[8],"ys") == 0 || strcmp(arg[9],"zs") == 0)
    error->all(FLERR,
               "Dump cfg arguments can not mix xs|ys|zs with xsu|ysu|zsu");
  else unwrapflag = 1;

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
          error->all(FLERR,"Invalid keyword in dump cfg command");
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
  if (auxname) {
    for (int i = 0; i < nfield-5; i++) delete [] auxname[i];
    delete [] auxname;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCFG::init_style()
{
  if (multifile == 0) 
    error->all(FLERR,"Dump cfg requires one snapshot per file");

  DumpCustom::init_style();

  // setup function ptrs

  if (buffer_flag == 1) write_choice = &DumpCFG::write_string;
  else write_choice = &DumpCFG::write_lines;
}

/* ---------------------------------------------------------------------- */

void DumpCFG::write_header(bigint n)
{
  // special handling for atom style peri
  //   use average volume of particles to scale particles to mimic C atoms
  //   scale box dimension to sc lattice for C with sigma = 1.44 Angstroms
  // special handling for unwrapped coordinates

  double scale;
  if (atom->peri_flag) {
    int nlocal = atom->nlocal;
    double vone = 0.0;
    for (int i = 0; i < nlocal; i++) vone += atom->vfrac[i];
    double vave;
    MPI_Allreduce(&vone,&vave,1,MPI_DOUBLE,MPI_SUM,world);
    if (atom->natoms) vave /= atom->natoms;
    if (vave > 0.0) scale = 1.44 / pow(vave,1.0/3.0);
  } else if (unwrapflag == 1) scale = UNWRAPEXPAND;
  else scale = 1.0;

  char str[64];
  sprintf(str,"Number of particles = %s\n",BIGINT_FORMAT);
  fprintf(fp,str,n);
  fprintf(fp,"A = %g Angstrom (basic length-scale)\n",scale);
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

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpCFG::convert_string(int n, double *mybuf)
{
  int i,j;

  int offset = 0;
  int m = 0;

  if (unwrapflag == 0) {
    for (i = 0; i < n; i++) {
      if (offset + size_one*ONEFIELD > maxsbuf) {
        if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
        maxsbuf += DELTA;
        memory->grow(sbuf,maxsbuf,"dump:sbuf");
      }

      for (j = 0; j < size_one; j++) {
        if (j == 0) {
	  offset += sprintf(&sbuf[offset],"%f \n",mybuf[m]);
        } else if (j == 1) {
	  offset += sprintf(&sbuf[offset],"%s \n",typenames[(int) mybuf[m]]);
        } else if (j >= 2) {
          if (vtype[j] == INT) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
          else if (vtype[j] == DOUBLE) 
            offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
          else if (vtype[j] == STRING) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],typenames[(int) mybuf[m]]);
          else if (vtype[j] == BIGINT) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],static_cast<bigint> (mybuf[m]));
        }
        m++;
      }
      offset += sprintf(&sbuf[offset],"\n");
    }

  } else if (unwrapflag == 1) {
    double unwrap_coord;
    for (i = 0; i < n; i++) {
      if (offset + size_one*ONEFIELD > maxsbuf) {
        if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
        maxsbuf += DELTA;
        memory->grow(sbuf,maxsbuf,"dump:sbuf");
      }

      for (j = 0; j < size_one; j++) {
        if (j == 0) {
	  offset += sprintf(&sbuf[offset],"%f \n",mybuf[m]);
        } else if (j == 1) {
	  offset += sprintf(&sbuf[offset],"%s \n",typenames[(int) mybuf[m]]);
        } else if (j >= 2 && j <= 4) {
          unwrap_coord = (mybuf[m] - 0.5)/UNWRAPEXPAND + 0.5;
          offset += sprintf(&sbuf[offset],vformat[j],unwrap_coord);
        } else if (j >= 5 ) {
          if (vtype[j] == INT) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],static_cast<int> (mybuf[m]));
          else if (vtype[j] == DOUBLE) 
            offset += sprintf(&sbuf[offset],vformat[j],mybuf[m]);
          else if (vtype[j] == STRING) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],typenames[(int) mybuf[m]]);
          else if (vtype[j] == BIGINT) 
            offset += 
              sprintf(&sbuf[offset],vformat[j],static_cast<bigint> (mybuf[m]));
        }
        m++;
      }
      offset += sprintf(&sbuf[offset],"\n");
    }
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpCFG::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCFG::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCFG::write_lines(int n, double *mybuf)
{
  int i,j,m;

  if (unwrapflag == 0) {
    m = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < size_one; j++) {
        if (j == 0) {
          fprintf(fp,"%f \n",mybuf[m]);
        } else if (j == 1) {
          fprintf(fp,"%s \n",typenames[(int) mybuf[m]]);
        } else if (j >= 2) {
          if (vtype[j] == INT) 
            fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
          else if (vtype[j] == DOUBLE) 
            fprintf(fp,vformat[j],mybuf[m]);
          else if (vtype[j] == STRING) 
            fprintf(fp,vformat[j],typenames[(int) mybuf[m]]);
          else if (vtype[j] == BIGINT) 
            fprintf(fp,vformat[j],static_cast<bigint> (mybuf[m]));
        }
        m++;
      }
      fprintf(fp,"\n");
    }
  } else if (unwrapflag == 1) {
    m = 0;
    double unwrap_coord;
    for (i = 0; i < n; i++) {
      for (j = 0; j < size_one; j++) {
        if (j == 0) {
          fprintf(fp,"%f \n",mybuf[m]);
        } else if (j == 1) {
          fprintf(fp,"%s \n",typenames[(int) mybuf[m]]);
        } else if (j >= 2 && j <= 4) {
          unwrap_coord = (mybuf[m] - 0.5)/UNWRAPEXPAND + 0.5;
          fprintf(fp,vformat[j],unwrap_coord);
        } else if (j >= 5 ) {
          if (vtype[j] == INT) 
            fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
          else if (vtype[j] == DOUBLE) 
            fprintf(fp,vformat[j],mybuf[m]);
          else if (vtype[j] == STRING) 
            fprintf(fp,vformat[j],typenames[(int) mybuf[m]]);
          else if (vtype[j] == BIGINT) 
            fprintf(fp,vformat[j],static_cast<bigint> (mybuf[m]));
        }
        m++;
      }
      fprintf(fp,"\n");
    }
  }
}
