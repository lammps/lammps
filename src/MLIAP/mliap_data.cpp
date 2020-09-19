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

#include "mliap_data.h"

#include "atom.h"
#include "memory.h"
#include "mliap_descriptor.h"
#include "mliap_model.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

MLIAPData::MLIAPData(LAMMPS *lmp, int gradgradflag_in, int *map_in,
                     class MLIAPModel* model_in,
                     class MLIAPDescriptor* descriptor_in,
                     class PairMLIAP* pairmliap_in) :
  Pointers(lmp), gradforce(nullptr), betas(nullptr), descriptors(nullptr), gamma(nullptr),
  gamma_row_index(nullptr), gamma_col_index(nullptr), egradient(nullptr),
  numneighs(nullptr), iatoms(nullptr), ielems(nullptr), jatoms(nullptr), jelems(nullptr),
  rij(nullptr), graddesc(nullptr), model(nullptr), descriptor(nullptr), list(nullptr)
{
  gradgradflag = gradgradflag_in;
  map = map_in;
  model = model_in;
  descriptor = descriptor_in;
  pairmliap = pairmliap_in;

  ndescriptors = descriptor->ndescriptors;
  nelements = descriptor->nelements;
  nparams = model->get_nparams();

  gamma_nnz = model->get_gamma_nnz(this);
  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nparams*nelements;
  zoffset = 2*yoffset;
  natoms_array = atom->natoms;
  size_array_rows = 1+ndims_force*natoms_array+ndims_virial;
  size_array_cols = nparams*nelements+1;
  size_gradforce = ndims_force*nparams*nelements;

  natoms_max = 0;
  natomneigh_max = 0;
  nneigh_max = 0;
  nmax = 0;
  natomgamma_max = 0;

}

/* ---------------------------------------------------------------------- */

MLIAPData::~MLIAPData()
{
  memory->destroy(betas);
  memory->destroy(descriptors);
  memory->destroy(gamma_row_index);
  memory->destroy(gamma_col_index);
  memory->destroy(gamma);
  memory->destroy(egradient);
  memory->destroy(gradforce);

  memory->destroy(iatoms);
  memory->destroy(ielems);
  memory->destroy(numneighs);
  memory->destroy(jatoms);
  memory->destroy(jelems);
  memory->destroy(rij);
  memory->destroy(graddesc);
}

/* ---------------------------------------------------------------------- */

void MLIAPData::init()
{
  memory->create(egradient,nelements*nparams,"MLIAPData:egradient");
}

/* ----------------------------------------------------------------------
   generate neighbor arrays
------------------------------------------------------------------------- */

void MLIAPData::generate_neighdata(NeighList* list_in, int eflag_in, int vflag_in)
{
  list = list_in;
  double **x = atom->x;
  int *type = atom->type;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // grow nmax gradforce array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(gradforce,nmax,size_gradforce,
                 "MLIAPData:gradforce");
 }

  // clear gradforce array

  int ntotal = atom->nlocal + atom->nghost;
  for (int i = 0; i < ntotal; i++)
    for (int j = 0; j < size_gradforce; j++) {
      gradforce[i][j] = 0.0;
    }

  // grow gamma arrays if necessary

  if (gradgradflag == 1) {
    const int natomgamma = list->inum;
    if (natomgamma_max < natomgamma) {
      memory->grow(gamma_row_index,natomgamma,gamma_nnz,"MLIAPData:gamma_row_index");
      memory->grow(gamma_col_index,natomgamma,gamma_nnz,"MLIAPData:gamma_col_index");
      memory->grow(gamma,natomgamma,gamma_nnz,"MLIAPData:gamma");
      natomgamma_max = natomgamma;
    }
  }

  // grow arrays if necessary

  natoms = list->inum;
  if (natoms_max < natoms) {
    memory->grow(betas,natoms,ndescriptors,"MLIAPData:betas");
    memory->grow(descriptors,natoms,ndescriptors,"MLIAPData:descriptors");
    natoms_max = natoms;
  }

  grow_neigharrays();

  int ij = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    const int i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];

      if (rsq < descriptor->cutsq[ielem][jelem]) {
        jatoms[ij] = j;
        jelems[ij] = jelem;
        rij[ij][0] = delx;
        rij[ij][1] = dely;
        rij[ij][2] = delz;
        ij++;
        ninside++;
      }
    }
    iatoms[ii] = i;
    ielems[ii] = ielem;
    numneighs[ii] = ninside;
  }

  eflag = eflag_in;
  vflag = vflag_in;
}

/* ----------------------------------------------------------------------
   grow neighbor arrays to handle all neighbors
------------------------------------------------------------------------- */

void MLIAPData::grow_neigharrays()
{

  // grow neighbor atom arrays if necessary

  const int natomneigh = list->inum;
  if (natomneigh_max < natomneigh) {
    memory->grow(iatoms,natomneigh,"MLIAPData:iatoms");
    memory->grow(ielems,natomneigh,"MLIAPData:ielems");
    memory->grow(numneighs,natomneigh,"MLIAPData:numneighs");
    natomneigh_max = natomneigh;
  }

  // grow neighbor arrays if necessary

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int iilast = list->inum-1;
  int ilast = list->ilist[iilast];
  int upperbound = firstneigh[ilast] - firstneigh[0] + numneigh[ilast];
  if (nneigh_max >= upperbound) return;

  double **x = atom->x;
  int *type = atom->type;

  int nneigh = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    const int i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    int ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const double delx = x[j][0] - xtmp;
      const double dely = x[j][1] - ytmp;
      const double delz = x[j][2] - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      const int jelem = map[jtype];
      if (rsq < descriptor->cutsq[ielem][jelem]) ninside++;
    }
    nneigh += ninside;
  }

  if (nneigh_max < nneigh) {
    memory->grow(jatoms,nneigh,"MLIAPData:jatoms");
    memory->grow(jelems,nneigh,"MLIAPData:jelems");
    memory->grow(rij,nneigh,3,"MLIAPData:rij");
    if (gradgradflag == 0)
      memory->grow(graddesc,nneigh,ndescriptors,3,"MLIAPData:graddesc");
    nneigh_max = nneigh;
  }
}

double MLIAPData::memory_usage()
{
  double bytes = 0.0;

  bytes += nelements*nparams*sizeof(double);     // egradient
  bytes += nmax*size_gradforce*sizeof(double);   // gradforce

  if (gradgradflag == 1) {
    bytes += natomgamma_max*
      gamma_nnz*sizeof(int);                     //gamma_row_index
    bytes += natomgamma_max*
      gamma_nnz*sizeof(int);                     // gamma_col_index
    bytes += natomgamma_max*
      gamma_nnz*sizeof(double);                  // gamma
  }

  bytes += natoms*ndescriptors*sizeof(int);      // betas
  bytes += natoms*ndescriptors*sizeof(int);      // descriptors

  bytes += natomneigh_max*sizeof(int);               // iatoms
  bytes += natomneigh_max*sizeof(int);               // ielems
  bytes += natomneigh_max*sizeof(int);               // numneighs

  bytes += nneigh_max*sizeof(int);                   // jatoms
  bytes += nneigh_max*sizeof(int);                   // jelems
  bytes += nneigh_max*3*sizeof(double);              // rij"

  if (gradgradflag == 0)
    bytes += nneigh_max*ndescriptors*3*sizeof(double);// graddesc

  return bytes;
}

