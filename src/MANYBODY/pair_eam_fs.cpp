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
   Contributing authors: Tim Lau (MIT)
------------------------------------------------------------------------- */

#include "pair_eam_fs.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairEAMFS::PairEAMFS(LAMMPS *lmp) : PairEAM(lmp)
{
  one_coeff = 1;
  manybody_flag = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read EAM Finnis-Sinclair file
------------------------------------------------------------------------- */

void PairEAMFS::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EAM Finnis-Sinclair file

  if (fs) {
    for (i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
  }
  fs = new Fs();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < fs->nelements; j++)
      if (strcmp(arg[i],fs->elements[j]) == 0) break;
    if (j < fs->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,fs->mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMFS::read_file(char *filename)
{
  Fs *file = fs;

  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
  int nwords = atom->count_words(line);
  if (nwords != file->nelements + 1)
    error->all(FLERR,"Incorrect element names in EAM potential file");

  char **words = new char*[file->nelements+1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

  file->elements = new char*[file->nelements];
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
    file->elements[i] = new char[n];
    strcpy(file->elements[i],words[i]);
  }
  delete [] words;

  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    nwords = sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&nwords,1,MPI_INT,0,world);
  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  if ((nwords != 5) || (file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
    error->all(FLERR,"Invalid EAM potential file");

  file->mass = new double[file->nelements];
  memory->create(file->frho,file->nelements,file->nrho+1,
                                              "pair:frho");
  memory->create(file->rhor,file->nelements,file->nelements,
                 file->nr+1,"pair:rhor");
  memory->create(file->z2r,file->nelements,file->nelements,
                 file->nr+1,"pair:z2r");

  int i,j,tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }
    MPI_Bcast(&file->mass[i],1,MPI_DOUBLE,0,world);

    if (me == 0) grab(fptr,file->nrho,&file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1],file->nrho,MPI_DOUBLE,0,world);

    for (j = 0; j < file->nelements; j++) {
      if (me == 0) grab(fptr,file->nr,&file->rhor[i][j][1]);
      MPI_Bcast(&file->rhor[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr,file->nr,&file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   copy read-in setfl potential to standard array format
------------------------------------------------------------------------- */

void PairEAMFS::file2array()
{
  int i,j,m,n;
  int ntypes = atom->ntypes;

  // set function params directly from fs file

  nrho = fs->nrho;
  nr = fs->nr;
  drho = fs->drho;
  dr = fs->dr;
  rhomax = (nrho-1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of fs elements + 1 for zero array

  nfrho = fs->nelements + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < fs->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = fs->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = square of # of fs elements

  nrhor = fs->nelements * fs->nelements;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // copy each element pair rhor to global rhor

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j < fs->nelements; j++) {
      for (m = 1; m <= nr; m++) rhor[n][m] = fs->rhor[i][j][m];
      n++;
    }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for fs files, there is a full NxN set of rhor arrays
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i] * fs->nelements + map[j];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of fs elements

  nz2r = fs->nelements * (fs->nelements+1) / 2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = fs->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}
