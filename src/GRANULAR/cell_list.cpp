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

#include "cell_list.h"
#include <assert.h>
#include "domain.h"
#include "error.h"
#include "comm.h"
#include <algorithm>
#include <numeric>
#include <mpi.h>

#define SMALL 1.0e-6

using namespace LAMMPS_NS;
using namespace std;

CellList::CellList(LAMMPS *lmp) : Pointers(lmp)
{
}

/* ----------------------------------------------------------------------
   insert element into cell list
------------------------------------------------------------------------- */

void CellList::insert(double * x, double r)
{
  int ibin = coord2bin(x);
  assert(ibin >= 0 && ibin <= nbins);

  int new_idx = elements.size();
  next.push_back(binhead[ibin]);
  elements.push_back(CellList::Element(x, r));

  binhead[ibin] = new_idx;
}

/* ----------------------------------------------------------------------
   remove all elements and reset cell list
------------------------------------------------------------------------- */

void CellList::clear()
{
  fill(binhead.begin(), binhead.end(), -1);
  next.clear();
  elements.clear();
}

/* ----------------------------------------------------------------------
   check if element at position x with radius r would overlap with any
   previously inserted element of cell list
------------------------------------------------------------------------- */

bool CellList::has_overlap(double *x, double r) const
{
  int ibin = coord2bin(x);
  if (ibin < 0 || ibin >= nbins) return false;

  for (int k = 0; k < NSTENCIL; ++k) {
    const int offset = stencil[k];
    for (int j = binhead[ibin+stencil[k]]; j >= 0; j = next[j]) {
        const CellList::Element & elem = elements[j];
        double dx = x[0] - elem.x[0];
        double dy = x[1] - elem.x[1];
        double dz = x[2] - elem.x[2];
        domain->minimum_image(dx, dy, dz);
        const double rsq = dx*dx + dy*dy + dz*dz;
        const double radsum = r + elem.r;
        if (rsq <= radsum*radsum) return true;
    }
  }

  return false;
}

/* ----------------------------------------------------------------------
   setup cell list for bounding box and binsize
------------------------------------------------------------------------- */

void CellList::setup(double * bboxlo, double * bboxhi, double binsize)
{
  double extent[3];

  // determine extent of bounding box
  extent[0] = bboxhi[0] - bboxlo[0];
  extent[1] = bboxhi[1] - bboxlo[1];
  extent[2] = bboxhi[2] - bboxlo[2];

  // save bounding box
  this->bboxlo[0] = bboxlo[0];
  this->bboxlo[1] = bboxlo[1];
  this->bboxlo[2] = bboxlo[2];
  this->bboxhi[0] = bboxhi[0];
  this->bboxhi[1] = bboxhi[1];
  this->bboxhi[2] = bboxhi[2];

  const double binsizeinv = 1.0/binsize;

  // test for too many bins in any dimension due to huge global domain

  assert(extent[0]*binsizeinv < MAXSMALLINT &&
         extent[1]*binsizeinv < MAXSMALLINT &&
         extent[2]*binsizeinv < MAXSMALLINT);

  // create actual bins

  nbinx = static_cast<int> (extent[0]*binsizeinv);
  nbiny = static_cast<int> (extent[1]*binsizeinv);
  nbinz = static_cast<int> (extent[2]*binsizeinv);

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly

  binsizex = extent[0]/nbinx;
  binsizey = extent[1]/nbiny;
  binsizez = extent[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  // mbinlo/hi = lowest and highest global bins
  // coord = lowest and highest values of coords
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  double * bsubboxlo = bboxlo;
  double * bsubboxhi = bboxhi;
  int mbinxhi,mbinyhi,mbinzhi;
  double coord;

  coord = bsubboxlo[0] - SMALL*extent[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*extent[0];
  mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*extent[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*extent[1];
  mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

  coord = bsubboxlo[2] - SMALL*extent[2];
  mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
  if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
  coord = bsubboxhi[2] + SMALL*extent[2];
  mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);

  // extend bins by 1 to insure stencil extent is included

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  mbinz = mbinzhi - mbinzlo + 1;

  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz) + 1;
  assert(bbin < MAXSMALLINT);

  nbins = bbin;
  binhead.resize(bbin);
  fill(binhead.begin(), binhead.end(), -1);

  // create 3D stencil for cell lookup
  int s = 0;

  for (int k = -1; k <= 1; ++k) {
    for (int j = -1; j <= 1; ++j) {
      for (int i = -1; i <= 1; ++i) {
        stencil[s++] = k*mbinx*mbiny + j*mbinx + i;
      }
    }
  }
}

int CellList::coord2bin(double *x) const
{
  int ix, iy, iz;

  if (x[0] >= bboxhi[0]) {
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  } else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = std::min(ix,nbinx-1);
  } else {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;
  }

  if (x[1] >= bboxhi[1]) {
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  } else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = std::min(iy,nbiny-1);
  } else {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;
  }

  if (x[2] >= bboxhi[2]) {
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  } else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = std::min(iz,nbinz-1);
  } else {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;
  }

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}


size_t CellList::count() const {
    return elements.size();
}

/* -------------------------------------------------------------------------- */

DistributedCellList::DistributedCellList(LAMMPS * lmp) : CellList(lmp) {
  recvcounts.resize(comm->nprocs);
  displs.resize(comm->nprocs);

  MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
  int blocklenghts[2]   = {3, 1};
  MPI_Aint offsets[2] = {0, 3*sizeof(double)};
  MPI_Type_create_struct(2, blocklenghts, offsets, types, &mpi_element_type);
  MPI_Type_commit(&mpi_element_type);
}

void DistributedCellList::allgather(INearList * local_nlist) {
  CellList * clist = dynamic_cast<CellList*>(local_nlist);

  if(!clist) {
    error->all(FLERR,"DistributedCellList::allgather requires pointer to CellList object!");
  }


  int ncount_local = clist->count();
  assert(nbins == clist->nbins);

  MPI_Allgather(&ncount_local, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, world);
  int ncount_global = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

  const int nprocs = comm->nprocs;

  binhead.resize(nprocs * nbins);
  next.resize(ncount_global);
  elements.resize(ncount_global);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++) {
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
  }

  // collect binheads from all processors
  MPI_Allgather(&clist->binhead[0], nbins, MPI_INT, &binhead[0], nbins, MPI_INT, world);

  // collect next array from all processors that have elements
  int * nextptr = NULL;
  if (ncount_local) nextptr = &clist->next[0];
  MPI_Allgatherv(nextptr, ncount_local, MPI_INT, &next[0], &recvcounts[0], &displs[0], MPI_INT, world);

  // collect element arrays from all processors that have elements
  CellList::Element * elementptr = NULL;
  if (ncount_local) elementptr = &clist->elements[0];
  MPI_Allgatherv(elementptr, ncount_local, mpi_element_type, &elements[0], &recvcounts[0], &displs[0], mpi_element_type, world);

  // at this point all processors have a full list of all elements in the global cell list,
  // however, their binheads still need to be merged

  // step 1: correct element indices, each process has been counting from 0.
  //         now that we have them in one array we can reuse the displs array to correct the indices

  for (int iproc = 1; iproc < nprocs; iproc++) {
      const int offset = displs[iproc];
      const int ibin_offset = iproc * nbins;

      // correct binhead index
      for(int ibin = 0; ibin < nbins; ++ibin) {
          int & head = binhead[ibin+ibin_offset];
          if (head >= 0) {
              head += offset;
          }
      }

      // correct next index
      for(int i = 0; i < recvcounts[iproc]; ++i) {
          int & n = next[offset + i];
          if (n >= 0) {
              n += offset;
          }
      }
  }

  // step 2: reduce binheads and next to single binning structure by merging lists
  for(int ibin = 0; ibin < nbins; ++ibin) {
      int iproc = 0;
      int j = binhead[ibin];
      int prev = -1;

      while(iproc < nprocs) {
          while(j >= 0) {
              prev = j;
              j = next[j];
          }

          iproc++;

          if(iproc < nprocs) {
              const int ibin_offset = iproc * nbins;
              j = binhead[ibin_offset + ibin];
              if(prev >= 0) {
                next[prev] = j;
              }
          }
      }
  }
}
