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
#include <unistd.h>
#include "write_coeff.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "comm.h"
#include "force.h"
#include "universe.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   called as write_coeff command in input script
------------------------------------------------------------------------- */

void WriteCoeff::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_coeff command before simulation box is defined");

  if (narg != 1) error->all(FLERR,"Illegal write_coeff command");

  int n = strlen(arg[0]) + 5;
  char *file = new char[n];

  strcpy(file,"tmp.");
  strcat(file,arg[0]);

  // initialize relevant styles
  force->init();

  if (comm->me == 0) {
    char str[256], coeff[256];
    FILE *one = fopen(file,"wb+");
    if (one == NULL) {
          sprintf(str,"Cannot open coeff file %s",file);
      error->one(FLERR,str);
    }

    if (force->pair && force->pair->writedata) {
      fprintf(one,"# pair_style %s\npair_coeff\n",force->pair_style);
      force->pair->write_data_all(one);
      fprintf(one,"end\n");
    }
    if (force->bond && force->bond->writedata) {
      fprintf(one,"# bond_style %s\nbond_coeff\n",force->bond_style);
      force->bond->write_data(one);
      fprintf(one,"end\n");
    }
    if (force->angle && force->angle->writedata) {
      fprintf(one,"# angle_style %s\nangle_coeff\n",force->angle_style);
      force->angle->write_data(one);
      fprintf(one,"end\n");
    }
    if (force->dihedral && force->dihedral->writedata) {
      fprintf(one,"# dihedral_style %s\ndihedral_coeff\n",
              force->dihedral_style);
      force->dihedral->write_data(one);
      fprintf(one,"end\n");
    }
    if (force->improper && force->improper->writedata) {
      fprintf(one,"# improper_style %s\nimproper_coeff\n",
              force->improper_style);
      force->improper->write_data(one);
      fprintf(one,"end\n");
    }
    rewind(one);

    FILE *two = fopen(file+4,"w");
    if (two == NULL) {
      sprintf(str,"Cannot open coeff file %s",file+4);
      error->one(FLERR,str);
    }
    fprintf(two,"# LAMMPS coeff file via write_coeff, version %s\n",
            universe->version);
    while(1) {
      if (fgets(str,256,one) == NULL) break;
      fputs(str,two);      // style
      fgets(str,256,one);  // coeff
      n = strlen(str);
      strcpy(coeff,str);
      coeff[n-1] = '\0';
      fgets(str,256,one);
      while (strcmp(str,"end\n") != 0) {
        fprintf(two,"%s %s",coeff,str);
        fgets(str,256,one);
      }
      fputc('\n',two);
    }
    fclose(one);
    fclose(two);
    unlink(file);
  }

  delete [] file;
}
