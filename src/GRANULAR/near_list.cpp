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
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "near_list.h"
#include <assert.h>
#include "domain.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include <mpi.h>

using namespace LAMMPS_NS;
using namespace std;

NearList::NearList(LAMMPS *lmp) : Pointers(lmp), elements(nullptr), ncount(0), nallocated(0)
{
}

NearList::~NearList()
{
  // free local memory
  if(elements) {
    memory->destroy(elements);
  }
}

/* ----------------------------------------------------------------------
   insert element into near list
------------------------------------------------------------------------- */

void NearList::insert(double * x, double r)
{
  elements[ncount][0] = x[0];
  elements[ncount][1] = x[1];
  elements[ncount][2] = x[2];
  elements[ncount][3] = r;
  ncount++;
}

void NearList::allocate(int nmax){
  memory->create(elements, nmax, 4, "nearlist:elements");
  nallocated = nmax;
}

/* ----------------------------------------------------------------------
   check if element at position x with radius r would overlap with any
   previously inserted element of cell list
------------------------------------------------------------------------- */

bool NearList::has_overlap(double *x, double r) const
{
  for (int i = 0; i < ncount; i++) {
    double dx = x[0] - elements[i][0];
    double dy = x[1] - elements[i][1];
    double dz = x[2] - elements[i][2];

    // use minimum_image() to account for PBC
    domain->minimum_image(dx, dy, dz);

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double radsum = r + elements[i][3];

    if (rsq <= radsum*radsum) {
      // found overlap
      return true;
    }
  }
  return false;
}

size_t NearList::count() const {
  return ncount;
}



DistributedNearList::DistributedNearList(LAMMPS * lmp) : NearList(lmp) {
  // allgather arrays
  recvcounts = new int[comm->nprocs];
  displs = new int[comm->nprocs];
}

DistributedNearList::~DistributedNearList() {
  delete [] recvcounts;
  delete [] displs;
}

void DistributedNearList::allgather(INearList * local_nlist) {
  NearList * nlist = dynamic_cast<NearList*>(local_nlist);

  if(!nlist) {
    error->all(FLERR,"DistributedNearList::allgather requires pointer to NearList object!");
  }

  assert(nallocated >= nlist->count());

  // setup for allgatherv
  ncount = nlist->count();
  int n = 4 * ncount;
  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);

  int ntotal = 0;
  MPI_Allreduce(&ncount, &ntotal, 1, MPI_INT, MPI_SUM, world);

  displs[0] = 0;

  for (int iproc = 1; iproc < comm->nprocs; iproc++) {
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
  }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = nullptr;
  if (ncount) ptr = nlist->elements[0];
  MPI_Allgatherv(ptr, 4 * ncount, MPI_DOUBLE, elements[0], recvcounts, displs, MPI_DOUBLE, world);

  ncount = ntotal;
}
