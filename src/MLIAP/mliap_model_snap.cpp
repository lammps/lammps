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

#include "mliap_model_snap.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModelSNAP::MLIAPModelSNAP(LAMMPS* lmp) : Pointers(lmp)
{
  nelements = 0;
  elements = NULL;
  coeffelem = NULL;
  map = NULL;

  beta_max = 0;
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

MLIAPModelSNAP::~MLIAPModelSNAP()
{
  if (nelements) {
    for (int i = 0; i < nelements; i++)
      delete[] elements[i];
    delete[] elements;
    memory->destroy(coeffelem);
  }

  if (allocated)
    memory->destroy(map);

}

/* ----------------------------------------------------------------------
   Calculate model gradients for each atom i.e. dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelSNAP::gradient(NeighList* list, double **bispectrum, double *atomenergy, double **beta, int eflag)
{
  int i;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = map[itype];
    double* coeffi = coeffelem[ielem];

    for (int icoeff = 0; icoeff < ncoeff; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    if (quadraticflag) {
      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        double bveci = bispectrum[ii][icoeff];
        beta[ii][icoeff] += coeffi[k]*bveci;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          double bvecj = bispectrum[ii][jcoeff];
          beta[ii][icoeff] += coeffi[k]*bvecj;
          beta[ii][jcoeff] += coeffi[k]*bveci;
          k++;
        }
      }
    }

    if (eflag) {

      // energy of atom I, sum over coeffs_k * Bi_k

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E = beta.B + 0.5*B^t.alpha.B

      // linear contributions

      for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        etmp += coeffi[icoeff+1]*bispectrum[ii][icoeff];

      // quadratic contributions

      if (quadraticflag) {
        int k = ncoeff+1;
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          double bveci = bispectrum[ii][icoeff];
          etmp += 0.5*coeffi[k++]*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            double bvecj = bispectrum[ii][jcoeff];
            etmp += coeffi[k++]*bveci*bvecj;
          }
        }
      }
      atomenergy[ii] = etmp;
    }


  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void MLIAPModelSNAP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(map,n+1,"mliap_model_snap:map");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void MLIAPModelSNAP::init(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  if (nelements) {
    for (int i = 0; i < nelements; i++)
      delete[] elements[i];
    delete[] elements;
    memory->destroy(coeffelem);
  }

  char* type1 = arg[0];
  char* type2 = arg[1];
  char* coefffilename = arg[2];
  char* paramfilename = arg[3];
  char** elemtypes = &arg[4];

  // insure I,J args are * *

  if (strcmp(type1,"*") != 0 || strcmp(type2,"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read snapcoeff and snapparam files

  read_files(coefffilename,paramfilename);

  if (!quadraticflag)
    ncoeff = ncoeffall - 1;
  else {

    // ncoeffall should be (ncoeff+2)*(ncoeff+1)/2
    // so, ncoeff = floor(sqrt(2*ncoeffall))-1

    ncoeff = sqrt(2*ncoeffall)-1;
    ncoeffq = (ncoeff*(ncoeff+1))/2;
    int ntmp = 1+ncoeff+ncoeffq;
    if (ntmp != ncoeffall) {
      printf("ncoeffall = %d ntmp = %d ncoeff = %d \n",ncoeffall,ntmp,ncoeff);
      error->all(FLERR,"Incorrect SNAP coeff file");
    }
  }

  // read args that map atom types to SNAP elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  for (int i = 1; i <= atom->ntypes; i++) {
    char* elemname = elemtypes[i-1];
    int jelem;
    for (jelem = 0; jelem < nelements; jelem++)
      if (strcmp(elemname,elements[jelem]) == 0)
        break;

    if (jelem < nelements)
      map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) map[i] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

}

/* ---------------------------------------------------------------------- */

void MLIAPModelSNAP::read_files(char *coefffilename, char *paramfilename)
{

  // open SNAP coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = force->open_potential(coefffilename);
    if (fpcoeff == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open SNAP coefficient file %s",coefffilename);
      error->one(FLERR,str);
    }
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
    nwords = atom->count_words(line);
  }
  if (nwords != 2)
    error->all(FLERR,"Incorrect format in SNAP coefficient file");

  // words = ptrs to all words in line
  // strip single and double quotes from words

  char* words[MAXWORD];
  int iword = 0;
  words[iword] = strtok(line,"' \t\n\r\f");
  iword = 1;
  words[iword] = strtok(NULL,"' \t\n\r\f");

  nelements = atoi(words[0]);
  ncoeffall = atoi(words[1]);

  // set up element lists

  elements = new char*[nelements];
  memory->create(coeffelem,nelements,ncoeffall,"mliap_snap_model:coeffelem");

  // Loop over nelements blocks in the SNAP coefficient file

  for (int ielem = 0; ielem < nelements; ielem++) {

    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpcoeff);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof)
      error->all(FLERR,"Incorrect format in SNAP coefficient file");
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    nwords = atom->count_words(line);
    if (nwords != 3)
      error->all(FLERR,"Incorrect format in SNAP coefficient file");

    iword = 0;
    words[iword] = strtok(line,"' \t\n\r\f");
    iword = 1;
    words[iword] = strtok(NULL,"' \t\n\r\f");
    iword = 2;
    words[iword] = strtok(NULL,"' \t\n\r\f");

    char* elemtmp = words[0];
    int n = strlen(elemtmp) + 1;
    elements[ielem] = new char[n];
    strcpy(elements[ielem],elemtmp);

    for (int icoeff = 0; icoeff < ncoeffall; icoeff++) {
      if (comm->me == 0) {
        ptr = fgets(line,MAXLINE,fpcoeff);
        if (ptr == NULL) {
          eof = 1;
          fclose(fpcoeff);
        } else n = strlen(line) + 1;
      }

      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof)
        error->all(FLERR,"Incorrect format in SNAP coefficient file");
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);

      nwords = atom->count_words(line);
      if (nwords != 1)
        error->all(FLERR,"Incorrect format in SNAP coefficient file");

      iword = 0;
      words[iword] = strtok(line,"' \t\n\r\f");

      coeffelem[ielem][icoeff] = atof(words[0]);

    }
  }

  if (comm->me == 0) fclose(fpcoeff);

  // set flags for required keywords

  twojmaxflag = 0;

  // Set defaults for optional keywords

  quadraticflag = 0;
  alloyflag = 0;
  wselfallflag = 0;
  chunksize = 2000;

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = force->open_potential(paramfilename);
    if (fpparam == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open SNAP parameter file %s",paramfilename);
      error->one(FLERR,str);
    }
  }

  eof = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpparam);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpparam);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    if (nwords != 2)
      error->all(FLERR,"Incorrect format in SNAP parameter file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    char* keywd = strtok(line,"' \t\n\r\f");
    char* keyval = strtok(NULL,"' \t\n\r\f");

    if (comm->me == 0) {
      if (screen) fprintf(screen,"SNAP keyword %s %s \n",keywd,keyval);
      if (logfile) fprintf(logfile,"SNAP keyword %s %s \n",keywd,keyval);
    }

    if (strcmp(keywd,"twojmax") == 0) {
      twojmax = atoi(keyval);
      twojmaxflag = 1;
    } else if (strcmp(keywd,"quadraticflag") == 0)
      quadraticflag = atoi(keyval);
    else if (strcmp(keywd,"alloyflag") == 0)
      alloyflag = atoi(keyval);
    else if (strcmp(keywd,"wselfallflag") == 0)
      wselfallflag = atoi(keyval);
    else if (strcmp(keywd,"chunksize") == 0)
      chunksize = atoi(keyval);
    else if (strcmp(keywd,"rcutfac") == 0 || 
               strcmp(keywd,"rfac0") == 0 || 
               strcmp(keywd,"rmin0") == 0 ||
               strcmp(keywd,"switchflag") == 0 ||
               strcmp(keywd,"bzeroflag") == 0)
      continue;
    else
      error->all(FLERR,"Incorrect SNAP parameter file");
  }

  if (twojmaxflag == 0)
    error->all(FLERR,"Incorrect SNAP parameter file");

}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPModelSNAP::memory_usage()
{
  double bytes = 0;

  int n = atom->ntypes+1;
  bytes += n*sizeof(int);        // map

  return bytes;
}

