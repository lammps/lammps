// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "mliap_descriptor_snap.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data.h"
#include "pair_mliap.h"
#include "sna.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::MLIAPDescriptorSNAP(LAMMPS *lmp, char *paramfilename):
  MLIAPDescriptor(lmp)
{
  radelem = nullptr;
  wjelem = nullptr;
  snaptr = nullptr;
  read_paramfile(paramfilename);

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag, nelements);

  ndescriptors = snaptr->ncoeff;
}

/* ---------------------------------------------------------------------- */

MLIAPDescriptorSNAP::~MLIAPDescriptorSNAP()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  delete snaptr;
}

/* ----------------------------------------------------------------------
   compute descriptors for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptors(class MLIAPData* data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    // insure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      snaptr->element[ninside] = jelem; // element index for chem snap
      ninside++;
      ij++;
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

    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->descriptors[ii][icoeff] = snaptr->blist[icoeff];
  }

}

/* ----------------------------------------------------------------------
   compute forces for each atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_forces(class MLIAPData* data)
{
  double fij[3];
  double **f = atom->f;

  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];

    // insure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      snaptr->element[ninside] = jelem; // element index for chem snap
      ninside++;
      ij++;
    }

    // compute Ui, Yi for atom I

    if (chemflag)
      snaptr->compute_ui(ninside, ielem);
    else
      snaptr->compute_ui(ninside, 0);

    // for neighbors of I within cutoff:
    // compute Fij = dEi/dRj = -dEi/dRi
    // add to Fi, subtract from Fj

    snaptr->compute_yi(data->betas[ii]);

    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      if (chemflag)
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

      if (data->vflag)
        data->pairmliap->v_tally(i,j,fij,snaptr->rij[jj]);
    }
  }

}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_force_gradients(class MLIAPData* data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];

    // insure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];

      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      snaptr->element[ninside] = jelem; // element index for chem snap
      ninside++;
      ij++;
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

      if (chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_dbidrj();

      // Accumulate gamma_lk*dB_k/dRi, -gamma_lk**dB_k/dRj

      for (int inz = 0; inz < data->gamma_nnz; inz++) {
        const int l = data->gamma_row_index[ii][inz];
        const int k = data->gamma_col_index[ii][inz];
        data->gradforce[i][l]         += data->gamma[ii][inz]*snaptr->dblist[k][0];
        data->gradforce[i][l+data->yoffset] += data->gamma[ii][inz]*snaptr->dblist[k][1];
        data->gradforce[i][l+data->zoffset] += data->gamma[ii][inz]*snaptr->dblist[k][2];
        data->gradforce[j][l]         -= data->gamma[ii][inz]*snaptr->dblist[k][0];
        data->gradforce[j][l+data->yoffset] -= data->gamma[ii][inz]*snaptr->dblist[k][1];
        data->gradforce[j][l+data->zoffset] -= data->gamma[ii][inz]*snaptr->dblist[k][2];
      }

    }
  }

}

/* ----------------------------------------------------------------------
   compute descriptor gradients for each neighbor atom
   ---------------------------------------------------------------------- */

void MLIAPDescriptorSNAP::compute_descriptor_gradients(class MLIAPData* data)
{
  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    // insure rij, inside, wj, and rcutij are of size jnum

    const int jnum = data->numneighs[ii];
    snaptr->grow_rij(jnum);

    int ij0 = ij;
    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = data->jatoms[ij];
      const int jelem = data->jelems[ij];
      const double *delr = data->rij[ij];

      snaptr->rij[ninside][0] = delr[0];
      snaptr->rij[ninside][1] = delr[1];
      snaptr->rij[ninside][2] = delr[2];
      snaptr->inside[ninside] = j;
      snaptr->wj[ninside] = wjelem[jelem];
      snaptr->rcutij[ninside] = sqrt(cutsq[ielem][jelem]);
      snaptr->element[ninside] = jelem; // element index for chem snap
      ninside++;
      ij++;
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

    ij = ij0;
    for (int jj = 0; jj < ninside; jj++) {
      if (chemflag)
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, snaptr->element[jj]);
      else
        snaptr->compute_duidrj(snaptr->rij[jj], snaptr->wj[jj],
                               snaptr->rcutij[jj],jj, 0);

      snaptr->compute_dbidrj();

      // Accumulate dB_k^i/dRi, dB_k^i/dRj

      for (int k = 0; k < data->ndescriptors; k++) {
        data->graddesc[ij][k][0] = snaptr->dblist[k][0];
        data->graddesc[ij][k][1] = snaptr->dblist[k][1];
        data->graddesc[ij][k][2] = snaptr->dblist[k][2];
      }
      ij++;
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

  for (int i = 0; i < nelements; i++) delete[] elements[i];
  delete[] elements;
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);

  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = utils::open_potential(paramfilename,lmp,nullptr);
    if (fpparam == nullptr)
      error->one(FLERR,"Cannot open SNAP parameter file {}: {}",
                                   paramfilename, utils::getsyserror());
  }

  char line[MAXLINE],*ptr;
  int eof = 0;
  int n,nwords;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpparam);
      if (ptr == nullptr) {
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
    char* keyval = strtok(nullptr,"' \t\n\r\f");

    if (comm->me == 0)
      utils::logmesg(lmp,"SNAP keyword {} {} \n", keywd, keyval);

    // check for keywords with one value per element

    if (strcmp(keywd,"elems") == 0 ||
        strcmp(keywd,"radelems") == 0 ||
        strcmp(keywd,"welems") == 0) {

      if (nelementsflag == 0 || nwords != nelements+1)
        error->all(FLERR,"Incorrect SNAP parameter file");

      if (strcmp(keywd,"elems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          elements[ielem] = utils::strdup(keyval);
          keyval = strtok(nullptr,"' \t\n\r\f");
        }
        elementsflag = 1;
      } else if (strcmp(keywd,"radelems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          radelem[ielem] = utils::numeric(FLERR,keyval,false,lmp);
          keyval = strtok(nullptr,"' \t\n\r\f");
        }
        radelemflag = 1;
      } else if (strcmp(keywd,"welems") == 0) {
        for (int ielem = 0; ielem < nelements; ielem++) {
          wjelem[ielem] = utils::numeric(FLERR,keyval,false,lmp);
          keyval = strtok(nullptr,"' \t\n\r\f");
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
    for (int jelem = ielem+1; jelem < nelements; jelem++) {
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
  double bytes = MLIAPDescriptor::memory_usage();
  bytes += snaptr->memory_usage();                      // SNA object

  return bytes;
}
