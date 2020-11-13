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

#include "nbin.h"
#include <cmath>
#include "neighbor.h"
#include "neigh_request.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NBin::NBin(LAMMPS *lmp) : Pointers(lmp)
{
  last_bin = -1;
  mbins = maxbin = maxatom = 0;
  binhead = nullptr;
  bins = nullptr;
  atom2bin = nullptr;

  nbinx_multi2 = nullptr; nbiny_multi2 = nullptr; nbinz_multi2 = nullptr;
  mbins_multi2 = nullptr;
  mbinx_multi2 = nullptr; mbiny_multi2 = nullptr, mbinz_multi2 = nullptr;
  mbinxlo_multi2 = nullptr;
  mbinylo_multi2 = nullptr;
  mbinzlo_multi2 = nullptr;
  binsizex_multi2 = nullptr; binsizey_multi2 = nullptr; binsizez_multi2 = nullptr;
  bininvx_multi2 = nullptr; bininvy_multi2 = nullptr; bininvz_multi2 = nullptr;
  binhead_multi2 = nullptr;
  bins_multi2 = nullptr;
  atom2bin_multi2 = nullptr;
  maxbins_multi2 = nullptr;

  maxtypes = 0;

  neighbor->last_setup_bins = -1;

  // geometry settings

  dimension = domain->dimension;
  triclinic = domain->triclinic;

  kokkos = 0;
}

/* ---------------------------------------------------------------------- */

NBin::~NBin()
{
  memory->destroy(binhead);
  memory->destroy(bins);
  memory->destroy(atom2bin);
  
  if (!bins_multi2) return;
  
  memory->destroy(nbinx_multi2);
  memory->destroy(nbiny_multi2);
  memory->destroy(nbinz_multi2);
  memory->destroy(mbins_multi2);
  memory->destroy(mbinx_multi2);
  memory->destroy(mbiny_multi2);
  memory->destroy(mbinz_multi2);
  memory->destroy(mbinxlo_multi2);
  memory->destroy(mbinylo_multi2);
  memory->destroy(mbinzlo_multi2);

  memory->destroy(binsizex_multi2);
  memory->destroy(binsizey_multi2);
  memory->destroy(binsizez_multi2);
  memory->destroy(bininvx_multi2);
  memory->destroy(bininvy_multi2);
  memory->destroy(bininvz_multi2);

  for (int n = 1; n <= maxtypes; n++) {
    memory->destroy(binhead_multi2[n]);
    memory->destroy(bins_multi2[n]);
    memory->destroy(atom2bin_multi2[n]);
  }
  delete [] binhead_multi2;
  delete [] bins_multi2;
  delete [] atom2bin_multi2;

  memory->destroy(maxbins_multi2);  
}

/* ---------------------------------------------------------------------- */

void NBin::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class
------------------------------------------------------------------------- */

void NBin::copy_neighbor_info()
{
  includegroup = neighbor->includegroup;
  cutneighmin = neighbor->cutneighmin;
  cutneighmax = neighbor->cutneighmax;
  binsizeflag = neighbor->binsizeflag;
  binsize_user = neighbor->binsize_user;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

  // overwrite Neighbor cutoff with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) cutneighmax = cutoff_custom;
}



/* ----------------------------------------------------------------------
   convert atom coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NBin::coord2bin(double *x)
{
  int ix,iy,iz;

  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}


/* ----------------------------------------------------------------------
   convert atom coords into local bin # for a particular type
------------------------------------------------------------------------- */

int NBin::coord2bin_multi2(double *x, int it)
{
  int ix,iy,iz;
  int ibin;

  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx_multi2[it]) + nbinx_multi2[it];
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_multi2[it]);
    ix = MIN(ix,nbinx_multi2[it]-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_multi2[it]) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy_multi2[it]) + nbiny_multi2[it];
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_multi2[it]);
    iy = MIN(iy,nbiny_multi2[it]-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_multi2[it]) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz_multi2[it]) + nbinz_multi2[it];
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_multi2[it]);
    iz = MIN(iz,nbinz_multi2[it]-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_multi2[it]) - 1;

  
  ibin = (iz-mbinzlo_multi2[it])*mbiny_multi2[it]*mbinx_multi2[it]
       + (iy-mbinylo_multi2[it])*mbinx_multi2[it]
       + (ix-mbinxlo_multi2[it]);
  return ibin;
}

