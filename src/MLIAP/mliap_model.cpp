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

#include "mliap_model.h"
#include "pair_mliap.h"
#include <cstring>
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "utils.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "fmt/format.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModel::MLIAPModel(LAMMPS* lmp, char* coefffilename) : Pointers(lmp)
{
  coeffelem = NULL;
  if (coefffilename) read_coeffs(coefffilename);
  else {
    nparams = 0;
    nelements = 0;
    ndescriptors = 0;
  }
  nonlinearflag = 0;
}

/* ---------------------------------------------------------------------- */

MLIAPModel::~MLIAPModel()
{
  memory->destroy(coeffelem);
}

/* ----------------------------------------------------------------------
   placeholder
------------------------------------------------------------------------- */

void MLIAPModel::init()
{
}

/* ----------------------------------------------------------------------
   set number of elements
   ---------------------------------------------------------------------- */

void MLIAPModel::set_nelements(int nelements_in)
{
  nelements = nelements_in;
}

/* ----------------------------------------------------------------------
   set number of descriptors
   ---------------------------------------------------------------------- */

void MLIAPModel::set_ndescriptors(int ndescriptors_in)
{
  ndescriptors = ndescriptors_in;
}

/* ---------------------------------------------------------------------- */

void MLIAPModel::read_coeffs(char *coefffilename)
{

  // open coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = force->open_potential(coefffilename);
    if (fpcoeff == NULL)
      error->one(FLERR,fmt::format("Cannot open MLIAPModel coeff file {}: {}",
                                   coefffilename,utils::getsyserror()));
  }

  char line[MAXLINE],*ptr;
  int eof = 0;

  int n;
  int nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpcoeff);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = utils::count_words(line);
  }
  if (nwords != 2)
    error->all(FLERR,"Incorrect format in MLIAPModel coefficient file");

  // words = ptrs to all words in line
  // strip single and double quotes from words

  char* words[MAXWORD];
  int iword = 0;
  words[iword] = strtok(line,"' \t\n\r\f");
  iword = 1;
  words[iword] = strtok(NULL,"' \t\n\r\f");

  nelements = atoi(words[0]);
  nparams = atoi(words[1]);

  // set up coeff lists

  memory->create(coeffelem,nelements,nparams,"mliap_snap_model:coeffelem");

  // Loop over nelements blocks in the coefficient file

  for (int ielem = 0; ielem < nelements; ielem++) {
    for (int icoeff = 0; icoeff < nparams; icoeff++) {
      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fpcoeff);
        if (ptr == NULL) {
          eof = 1;
          fclose(fpcoeff);
        } else n = strlen(line) + 1;
      }

      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof)
        error->all(FLERR,"Incorrect format in  coefficient file");
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);

      nwords = utils::trim_and_count_words(line);
      if (nwords != 1)
        error->all(FLERR,"Incorrect format in  coefficient file");

      iword = 0;
      words[iword] = strtok(line,"' \t\n\r\f");

      coeffelem[ielem][icoeff] = atof(words[0]);

    }
  }

  if (comm->me == 0) fclose(fpcoeff);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPModel::memory_usage()
{
  double bytes = 0;

  bytes += nelements*nparams*sizeof(double);  // coeffelem
  return bytes;
}

