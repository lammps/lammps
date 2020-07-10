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
#include <cstdlib>
#include "mliap.h"
#include "mliap_model_linear.h"
#include "mliap_model_quadratic.h"
#include "mliap_descriptor_snap.h"
#include "compute_mliap.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

MLIAP::MLIAP(LAMMPS *lmp, int ndescriptors_in, int nparams_in, int nelements_in, 
             int gradgradflag_in, int *map_in, class MLIAPModel*, class MLIAPDescriptor*) :
  Pointers(lmp),
  list(NULL), 
  gradforce(NULL), 
  descriptors(NULL), gamma_row_index(NULL), gamma_col_index(NULL),
  gamma(NULL), egradient(NULL), model(NULL), descriptor(NULL),
  iatommliap(NULL), ielemmliap(NULL), numneighmliap(NULL),  
  jatommliap(NULL), jelemmliap(NULL), graddesc(NULL)
{
  ndescriptors = ndescriptors_in;
  nparams = nparams_in;
  nelements = nelements_in;
  gradgradflag = gradgradflag_in;
  map = map_in;

  gamma_nnz = model->get_gamma_nnz();
  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nparams*nelements;
  zoffset = 2*yoffset;
  natoms = atom->natoms;

  size_gradforce = ndims_force*nparams*nelements;

  natomdesc_max = 0;
  natomneigh_max = 0;
  nneigh_max = 0;

}

/* ---------------------------------------------------------------------- */

MLIAP::~MLIAP()
{
  memory->destroy(descriptors);
  memory->destroy(gamma_row_index);
  memory->destroy(gamma_col_index);
  memory->destroy(gamma);
  memory->destroy(egradient);
  memory->destroy(graddesc);

  memory->destroy(iatommliap);
  memory->destroy(ielemmliap);
  memory->destroy(numneighmliap);
  memory->destroy(jatommliap);
  memory->destroy(jelemmliap);
}

/* ---------------------------------------------------------------------- */

void MLIAP::init()
{
  memory->create(egradient,nelements*nparams,"ComputeMLIAP:egradient");
}

/* ----------------------------------------------------------------------
   generate neighbor arrays
------------------------------------------------------------------------- */

void MLIAP::generate_neigharrays(NeighList* list_in)
{
  list = list_in;
  double **x = atom->x;
  int *type = atom->type;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // grow arrays if necessary

  natomdesc = list->inum;
  if (natomdesc_max < natomdesc) {
    memory->grow(descriptors,natomdesc,ndescriptors,"ComputeMLIAP:descriptors");
    natomdesc_max = natomdesc;
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
        jatommliap[ij] = j;
        jelemmliap[ij] = jelem;
        ij++;
        ninside++;
      }
    }
    iatommliap[ii] = i;
    ielemmliap[ii] = ielem;
    numneighmliap[ii] = ninside;
  }
}

/* ----------------------------------------------------------------------
   grow neighbor arrays to handle all neighbors
------------------------------------------------------------------------- */

void MLIAP::grow_neigharrays()
{

  // grow neighbor atom arrays if necessary
    
  const int natomneigh = list->inum;
  if (natomneigh_max < natomneigh) {
    memory->grow(iatommliap,natomneigh,"ComputeMLIAP:iatommliap");
    memory->grow(ielemmliap,natomneigh,"ComputeMLIAP:ielemmliap");
    memory->grow(numneighmliap,natomneigh,"ComputeMLIAP:numneighmliap");
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
    memory->grow(jatommliap,nneigh,"ComputeMLIAP:jatommliap");
    memory->grow(jelemmliap,nneigh,"ComputeMLIAP:jelemmliap");
    memory->grow(graddesc,nneigh,ndescriptors,3,"ComputeMLIAP:graddesc");
    nneigh_max = nneigh;
  }
}

