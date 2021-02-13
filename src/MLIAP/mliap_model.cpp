/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "mliap_model.h"

#include "comm.h"
#include "error.h"
#include "memory.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModel::MLIAPModel(LAMMPS *lmp, char *) : Pointers(lmp)
{
  nparams = 0;
  nelements = 0;
  ndescriptors = 0;
  nlayers = 0;
  nonlinearflag = 0;
  nnflag = 0;
}

/* ---------------------------------------------------------------------- */

MLIAPModel::~MLIAPModel()
{
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


MLIAPModelSimple::MLIAPModelSimple(LAMMPS *lmp, char *coefffilename) : MLIAPModel(lmp, coefffilename)
{
  coeffelem = nullptr;
  nnodes = nullptr;
  activation = nullptr;
  scale = nullptr;
  if (coefffilename) read_coeffs(coefffilename);
}

/* ---------------------------------------------------------------------- */

MLIAPModelSimple::~MLIAPModelSimple()
{
  memory->destroy(coeffelem);
  if (nnflag) {
    memory->destroy(nnodes);
    memory->destroy(activation);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void MLIAPModelSimple::read_coeffs(char *coefffilename)
{

  // open coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = utils::open_potential(coefffilename,lmp,nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR,fmt::format("Cannot open MLIAPModel coeff file {}: {}",
                                   coefffilename,utils::getsyserror()));
  }

  char line[MAXLINE],*ptr, *tstr;
  int eof = 0;

  int n;
  int nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == nullptr) {
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
    if ((ptr = strstr(line,"nn"))) {
      *ptr = '\0';
      nnflag = 1;
    }
    nwords = utils::count_words(line);
  }
  if (nwords != 2)
    error->all(FLERR,"Incorrect format in MLIAPModel coefficient file");

  // words = ptrs to all words in line
  // strip single and double quotes from words

  try {
    ValueTokenizer coeffs(line);
    nelements = coeffs.next_int();
    nparams = coeffs.next_int();
  } catch (TokenizerException &e) {
    error->all(FLERR,fmt::format("Incorrect format in MLIAPModel coefficient "
                                 "file: {}",e.what()));
  }

  // set up coeff lists

  memory->create(coeffelem,nelements,nparams,"mliap_snap_model:coeffelem");

  if (nnflag == 0) {
      
  // Loop over nelements blocks in the coefficient file

    for (int ielem = 0; ielem < nelements; ielem++) {
      for (int icoeff = 0; icoeff < nparams; icoeff++) {
        if (comm->me == 0) {
          ptr = fgets(line,MAXLINE,fpcoeff);
          if (ptr == nullptr) {
            eof = 1;
            fclose(fpcoeff);
          } else n = strlen(line) + 1;
        }

        MPI_Bcast(&eof,1,MPI_INT,0,world);
        if (eof)
          error->all(FLERR,"Incorrect format in MLIAPModel coefficient file");
        MPI_Bcast(&n,1,MPI_INT,0,world);
        MPI_Bcast(line,n,MPI_CHAR,0,world);

        try {
          ValueTokenizer coeffs(utils::trim_comment(line));
          if (coeffs.count() != 1)
            throw TokenizerException("Wrong number of items","");
          coeffelem[ielem][icoeff] = coeffs.next_double();
        } catch (TokenizerException &e) {
          error->all(FLERR,fmt::format("Incorrect format in MLIAPModel "
                                     "coefficient file: {}",e.what()));
        }
      }
    }
    if (comm->me == 0) fclose(fpcoeff);
  }
    
  // set up the NN parameters
        
  else {
    int stats = 0;
    int ielem = 0;
    int l = 0;
          
    while (1) {
      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fpcoeff);
        if (ptr == nullptr) {
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
      nwords = utils::trim_and_count_words(line);
      if (nwords == 0) continue;
          
      if (stats == 0) { // Header NET
        tstr = strtok(line,"' \t\n\r\f");
        if (strncmp(tstr, "NET", 3) != 0) error->all(FLERR,"Incorrect format in NET coefficient file");
              
        ndescriptors = atoi(strtok(nullptr,"' \t\n\r\f"));
        nlayers = atoi(strtok(nullptr,"' \t\n\r\f"));
            
        memory->create(activation,nlayers,"mliap_model:activation");
        memory->create(nnodes,nlayers,"mliap_model:nnodes");
        memory->create(scale,nelements,2,ndescriptors,"mliap_model:scale");

        for (int ilayer = 0; ilayer < nlayers; ilayer++) {
          tstr = strtok(NULL,"' \t\n\r\f");
          nnodes[ilayer] = atoi(strtok(NULL,"' \t\n\r\f"));

          if (strncmp(tstr, "linear", 6) == 0) activation[ilayer] = 0;
          else if (strncmp(tstr, "sigmoid", 7) == 0) activation[ilayer] = 1;
          else if (strncmp(tstr, "tanh", 4) == 0) activation[ilayer] = 2;
          else if (strncmp(tstr, "relu", 4) == 0) activation[ilayer] = 3;
          else activation[ilayer] = 4;
        }
            
        stats = 1;
              
      } else if (stats == 1) {
          scale[ielem][0][l] = atof(strtok(line,"' \t\n\r\f"));
          for (int icoeff = 1; icoeff < nwords; icoeff++) {
            scale[ielem][0][l+icoeff] = atof(strtok(nullptr,"' \t\n\r\f"));
          }
          l += nwords;
          if (l == ndescriptors) {
            stats = 2;
            l = 0;
          }
              
      } else if (stats == 2) {
          scale[ielem][1][l] = atof(strtok(line,"' \t\n\r\f"));
          for (int icoeff = 1; icoeff < nwords; icoeff++) {
            scale[ielem][1][l+icoeff] = atof(strtok(nullptr,"' \t\n\r\f"));
          }
          l += nwords;
          if (l == ndescriptors) {
            stats = 3;
            l = 0;
          }
          
      // set up coeff lists
          
      } else if (stats == 3) {
     //         if (nwords > 30) error->all(FLERR,"Wrong number of items per line, max 30");
            
          coeffelem[ielem][l] = atof(strtok(line,"' \t\n\r\f"));
          for (int icoeff = 1; icoeff < nwords; icoeff++) {
            coeffelem[ielem][l+icoeff] = atof(strtok(nullptr,"' \t\n\r\f"));
          }
          l += nwords;
          if (l == nparams) {
            stats = 0;
            ielem++;
          }
      }
    }
    if (comm->me == 0) fclose(fpcoeff);
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPModelSimple::memory_usage()
{
  double bytes = 0;

  bytes += (double)nelements*nparams*sizeof(double);           // coeffelem
  if (nnflag) {
    bytes += (double)nelements*2*ndescriptors*sizeof(double);  // scale
    bytes += (int)nlayers*sizeof(int);                         // nnodes
    bytes += (int)nlayers*sizeof(int);                         // activation
  }
  return bytes;
}

