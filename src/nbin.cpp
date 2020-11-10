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

  nbinx_tiered = nullptr; nbiny_tiered = nullptr; nbinz_tiered = nullptr;
  mbins_tiered = nullptr;
  mbinx_tiered = nullptr; mbiny_tiered = nullptr, mbinz_tiered = nullptr;
  mbinxlo_tiered = nullptr;
  mbinylo_tiered = nullptr;
  mbinzlo_tiered = nullptr;
  binsizex_tiered = nullptr; binsizey_tiered = nullptr; binsizez_tiered = nullptr;
  bininvx_tiered = nullptr; bininvy_tiered = nullptr; bininvz_tiered = nullptr;
  binhead_tiered = nullptr;
  bins_tiered = nullptr;
  atom2bin_tiered = nullptr;
  maxbins_tiered = nullptr;

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
  
  if (!bins_tiered) return;
  
  memory->destroy(nbinx_tiered);
  memory->destroy(nbiny_tiered);
  memory->destroy(nbinz_tiered);
  memory->destroy(mbins_tiered);
  memory->destroy(mbinx_tiered);
  memory->destroy(mbiny_tiered);
  memory->destroy(mbinz_tiered);
  memory->destroy(mbinxlo_tiered);
  memory->destroy(mbinylo_tiered);
  memory->destroy(mbinzlo_tiered);

  memory->destroy(binsizex_tiered);
  memory->destroy(binsizey_tiered);
  memory->destroy(binsizez_tiered);
  memory->destroy(bininvx_tiered);
  memory->destroy(bininvy_tiered);
  memory->destroy(bininvz_tiered);

  for (int n = 1; n <= maxtypes; n++) {
    memory->destroy(binhead_tiered[n]);
    memory->destroy(bins_tiered[n]);
    memory->destroy(atom2bin_tiered[n]);
  }
  delete [] binhead_tiered;
  delete [] bins_tiered;
  delete [] atom2bin_tiered;

  memory->destroy(maxbins_tiered);  
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
