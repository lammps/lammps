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
