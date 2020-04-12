/* -*- c++ -*- ----------------------------------------------------------
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
#include <string>
#include <vector>
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "pair_species_coeff.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

void PairSpeciesCoeff::command(int narg, char **arg)
{
  if (narg != 1 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair species coefficients");

  // following the template of pair_comb read_file
  FILE *fp;
  fp = force->open_potential(arg[0]);
  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open potential file %s",arg[0]);
    error->one(FLERR,str);
  }

  std::vector<char *> species;
  for (int i = 0; i < atom->ntypes; ++i)
    species.push_back(arg[i+1]);

  char line[MAXLINE],*ptr;
  int n, eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    char *species1, *species2, *the_rest;
    ptr = line;
    species1 = strtok(ptr," \t");
    species2 = strtok(NULL," \t");
    the_rest = strtok(NULL,"\n");

    for (int type_a = 0; type_a < atom->ntypes; ++type_a) {
      for (int type_b = type_a; type_b < atom->ntypes; ++type_b) {
	if(((strcmp(species[type_a], species1) == 0) &&
            (strcmp(species[type_b], species2) == 0))
           ||
	   ((strcmp(species[type_b], species1) == 0) &&
            (strcmp(species[type_a], species2) == 0))
          ) {
          char pair_command[MAXLINE];
	  sprintf(pair_command, "pair_coeff %i %i %s", type_a+1, type_b+1,
                  the_rest);
	  input->one(pair_command);
	}
      }
    }
  }
  fclose(fp);
}
