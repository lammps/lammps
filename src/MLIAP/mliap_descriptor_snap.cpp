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

#include "mliap_descriptor_snap.h"
#include "pair_mliap.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "utils.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "sna.h"
#include "memory.h"
#include "error.h"
#include "fmt/format.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::MLIAPDescriptorSNAP(LAMMPS *lmp, char *paramfilename):
  MLIAPDescriptor(lmp)
{
  nelements = 0;
  elements = NULL;
  radelem = NULL;
  wjelem = NULL;
  snaptr = NULL;
  read_paramfile(paramfilename);

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag, nelements);

  ndescriptors = snaptr->ncoeff;

}

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::~MLIAPDescriptorSNAP()
{

  if (nelements) {
    for (int i = 0; i < nelements; i++)
      delete[] elements[i];
    delete[] elements;
    memory->destroy(radelem);
    memory->destroy(wjelem);
    memory->destroy(cutsq);
  }

  delete snaptr;

}

/* ----------------------------------------------------------------------
   compute descriptors for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptors(int* map, NeighList* list, double **descriptors)
{
  int i,j,jnum,ninside;
  double delx,dely,delz,rsq;
  int *jlist;

  double **x = atom->x;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];

    jlist = list->firstneigh[i];
    jnum = list->numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];

      if (rsq < cutsq[ielem][jelem]) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
	snaptr->wj[ninside] = wjelem[jelem];
	snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
	snaptr->element[ninside] = jelem; // element index for chem snap
        ninside++;
      }
    }

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);
    snaptr->compute_zi();

    if (chemflag)
      snaptr->compute_bi(ielem);
    else
      snaptr->compute_bi(0);

    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      descriptors[ii][icoeff] = snaptr->blist[icoeff];
  }

}

/* ----------------------------------------------------------------------
   compute forces for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_forces(PairMLIAP* pairmliap, NeighList* list, double **beta, int vflag)
{
  int i,j,jnum,ninside;
  double delx,dely,delz,rsq;
  double fij[3];
  int *jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = pairmliap->map[itype];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      int jelem = pairmliap->map[jtype];

      if (rsq < cutsq[ielem][jelem]) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
        snaptr->wj[ninside] = wjelem[jelem];
        snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
        snaptr->element[ninside] = jelem; // element index for chem snap
        ninside++;
      }
    }

    // compute Ui, Yi for atom I

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    // for neighbors of I within cutoff:
    // compute Fij = dEi/dRj = -dEi/dRi
    // add to Fi, subtract from Fj

    snaptr->compute_yi(beta[ii]);

    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      if(chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_deidrj(fij);

      f[i][0] += fij[0];
      f[i][1] += fij[1];
      f[i][2] += fij[2];
      f[j][0] -= fij[0];
      f[j][1] -= fij[1];
      f[j][2] -= fij[2];

      // add in global and per-atom virial contributions
      // this is optional and has no effect on force calculation

      if (vflag)
        pairmliap->v_tally(i,j,
                     fij[0],fij[1],fij[2],
                     -snaptr->rij[jj][0],-snaptr->rij[jj][1],
                     -snaptr->rij[jj][2]);

    }
  }

}

/* ----------------------------------------------------------------------
   compute force gradient for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_gradients(int *map, NeighList* list, 
                                         int gamma_nnz, int **gamma_row_index, 
                                         int **gamma_col_index, double **gamma, double **gradforce,
                                         int yoffset, int zoffset)
{
  int i,j,jnum,ninside;
  double delx,dely,delz,evdwl,rsq;
  double fij[3];
  int *jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];

      if (rsq < cutsq[ielem][jelem]) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
	snaptr->wj[ninside] = wjelem[jelem];
	snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
        snaptr->element[ninside] = jelem; // element index for chem snap
        ninside++;
      }
    }

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    snaptr->compute_zi();
    if (chemflag)
      snaptr->compute_bi(ielem);
    else
      snaptr->compute_bi(0);

    for (int jj = 0; jj < ninside; jj++) {
      const int j = snaptr->inside[jj];

      if(chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_dbidrj();
      
      // Accumulate gamma_lk*dB_k/dRi, -gamma_lk**dB_k/dRj
      
      for (int inz = 0; inz < gamma_nnz; inz++) {
        const int l = gamma_row_index[ii][inz];
        const int k = gamma_col_index[ii][inz];
        gradforce[i][l]         += gamma[ii][inz]*snaptr->dblist[k][0];
        gradforce[i][l+yoffset] += gamma[ii][inz]*snaptr->dblist[k][1];
        gradforce[i][l+zoffset] += gamma[ii][inz]*snaptr->dblist[k][2];
        gradforce[j][l]         -= gamma[ii][inz]*snaptr->dblist[k][0];
        gradforce[j][l+yoffset] -= gamma[ii][inz]*snaptr->dblist[k][1];
        gradforce[j][l+zoffset] -= gamma[ii][inz]*snaptr->dblist[k][2];
      }
      
    }
  }

}

/* ----------------------------------------------------------------------
   compute descriptor gradients for each neighbor atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptor_gradients(int *map, NeighList* list, 
                                         int gamma_nnz, int **gamma_row_index, 
                                         int **gamma_col_index, double **gamma, double **graddesc,
                                         int yoffset, int zoffset)
{
  int i,j,jnum,ninside;
  double delx,dely,delz,evdwl,rsq;
  double fij[3];
  int *jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum

    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];

      if (rsq < cutsq[ielem][jelem]) {
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
	snaptr->wj[ninside] = wjelem[jelem];
	snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
        snaptr->element[ninside] = jelem; // element index for chem snap
        ninside++;
      }
    }

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    snaptr->compute_zi();
    if (chemflag)
      snaptr->compute_bi(ielem);
    else
      snaptr->compute_bi(0);

    for (int jj = 0; jj < ninside; jj++) {
      const int j = snaptr->inside[jj];

      if(chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_dbidrj();

      // Accumulate dB_k^i/dRi, dB_k^i/dRj

      for (int k = 0; k < ndescriptors; k++) {
        graddesc[i][k] = snaptr->dblist[k][0];
        graddesc[i][k] = snaptr->dblist[k][1];
        graddesc[i][k] = snaptr->dblist[k][2];
        graddesc[j][k] = -snaptr->dblist[k][0];
        graddesc[j][k] = -snaptr->dblist[k][1];
        graddesc[j][k] = -snaptr->dblist[k][2];
      } 
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::init()
{
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::read_paramfile(char *paramfilename)
{

  // set flags for required keywords

  int rcutfacflag = 0;
  int twojmaxflag = 0;
  int nelementsflag = 0;
  int elementsflag = 0;
  int radelemflag = 0;
  int wjelemflag = 0;

  // Set defaults for optional keywords

  rfac0 = 0.99363;
  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = force->open_potential(paramfilename);
    if (fpparam == NULL)
      error->one(FLERR,fmt::format("Cannot open SNAP parameter file {}: {}",
                                   paramfilename, utils::getsyserror()));
  }

  char line[MAXLINE],*ptr;
  int eof = 0;
  int n,nwords;

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
    nwords = utils::count_words(line);
    if (nwords == 0) continue;

    // words = ptrs to all words in line
    // strip single and double quotes from words

    char* keywd = strtok(line,"' \t\n\r\f");
    char* keyval = strtok(NULL,"' \t\n\r\f");

    if (comm->me == 0) {
      utils::logmesg(lmp, fmt::format("SNAP keyword {} {} \n", keywd, keyval));
    }

    // check for keywords with one value per element

    if (strcmp(keywd,"elems") == 0 ||
        strcmp(keywd,"radelems") == 0 ||
        strcmp(keywd,"welems") == 0) {

      if (nelementsflag == 0 || nwords != nelements+1)
        error->all(FLERR,"Incorrect SNAP parameter file");

      if (strcmp(keywd,"elems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          char* elemtmp = keyval;
          int n = strlen(elemtmp) + 1;
          elements[ielem] = new char[n];
          strcpy(elements[ielem],elemtmp);
          keyval = strtok(NULL,"' \t\n\r\f");
        }
        elementsflag = 1;
      } else if (strcmp(keywd,"radelems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          radelem[ielem] = atof(keyval);
          keyval = strtok(NULL,"' \t\n\r\f");
        }
        radelemflag = 1;
      } else if (strcmp(keywd,"welems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          wjelem[ielem] = atof(keyval);
          keyval = strtok(NULL,"' \t\n\r\f");
        }
        wjelemflag = 1;
      }

    } else {

    // all other keywords take one value

      if (nwords != 2)
        error->all(FLERR,"Incorrect SNAP parameter file");

      if (strcmp(keywd,"nelems") == 0) {
        nelements = atoi(keyval);
        elements = new char*[nelements];
        memory->create(radelem,nelements,"mliap_snap_descriptor:radelem");
        memory->create(wjelem,nelements,"mliap_snap_descriptor:wjelem");
        nelementsflag = 1;
      } else if (strcmp(keywd,"rcutfac") == 0) {
        rcutfac = atof(keyval);
        rcutfacflag = 1;
      } else if (strcmp(keywd,"twojmax") == 0) {
        twojmax = atoi(keyval);
        twojmaxflag = 1;
      } else if (strcmp(keywd,"rfac0") == 0)
        rfac0 = atof(keyval);
      else if (strcmp(keywd,"rmin0") == 0)
        rmin0 = atof(keyval);
      else if (strcmp(keywd,"switchflag") == 0)
        switchflag = atoi(keyval);
      else if (strcmp(keywd,"bzeroflag") == 0)
        bzeroflag = atoi(keyval);
      else if (strcmp(keywd,"chemflag") == 0)
        chemflag = atoi(keyval);
      else if (strcmp(keywd,"bnormflag") == 0)
        bnormflag = atoi(keyval);
      else if (strcmp(keywd,"wselfallflag") == 0)
        wselfallflag = atoi(keyval);
      else
        error->all(FLERR,"Incorrect SNAP parameter file");

    }
  }

  if (!rcutfacflag || !twojmaxflag || !nelementsflag ||
      !elementsflag || !radelemflag || !wjelemflag)
    error->all(FLERR,"Incorrect SNAP parameter file");

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,nelements,nelements,"mliap/descriptor/snap:cutsq");
  for (int ielem = 0; ielem < nelements; ielem++) {
    cut = 2.0*radelem[ielem]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[ielem][ielem] = cut*cut;
    for(int jelem = ielem+1; jelem < nelements; jelem++) {
      cut = (radelem[ielem]+radelem[jelem])*rcutfac;
      cutsq[ielem][jelem] = cutsq[jelem][ielem] = cut*cut;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPDescriptorSNAP::memory_usage()
{
  double bytes = 0;

  bytes += snaptr->memory_usage(); // SNA object

  return bytes;
}

