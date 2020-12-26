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

#include "nbin_multi.h"
#include "neighbor.h"
#include "atom.h"
#include "group.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100

/* ---------------------------------------------------------------------- */

NBinMulti::NBinMulti(LAMMPS *lmp) : NBin(lmp) {}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinMulti::bin_atoms_setup(int nall) 
{
  // binhead_multi[n] = per-bin vector mbins in length mbins_multi[n]
  
  for (int n = 0; n < maxgroups; n++) {
    if (mbins_multi[n] > maxbins_multi[n]) {
      maxbins_multi[n] = mbins_multi[n];
      memory->destroy(binhead_multi[n]);
      memory->create(binhead_multi[n], mbins_multi[n], "neigh:mbins_multi");
    }
  }

  // bins_multi[n] and atom2bin_multi[n] = per-atom vectors
  // for both local and ghost atoms
  
  // TODO: maxatom should be maxatom per group

  if (nall > maxatom) {
    maxatom = nall;
    for (int n = 0; n < maxgroups; n++) {
      memory->destroy(bins_multi[n]);
      memory->destroy(atom2bin_multi[n]);
      memory->create(bins_multi[n], maxatom, "neigh:bin_multi");
      memory->create(atom2bin_multi[n], maxatom, "neigh:atom2bin_multi");
    }
  }
}

/* ---------------------------------------------------------------------
   Identify index of group with smallest cutoff
------------------------------------------------------------------------ */

int NBinMulti::igroup_min() 
{
  int imin = 0;
  for (int n = 0; n < maxgroups; n++)
    if (cutmultisq[n][n] < cutmultisq[imin][imin])  imin = n;

  return imin;
}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   ---------------------------------------------------------------------- */

void NBinMulti::setup_bins(int style)
{
  int n;
  int igroupmin;

  // Initialize arrays
  
  if (n_multi_groups > maxgroups) {

    // Clear any/all memory for existing types

    for (n = 0; n < maxgroups; n++) {
      memory->destroy(atom2bin_multi[n]);
      memory->destroy(binhead_multi[n]);
      memory->destroy(bins_multi[n]);
    }
    delete [] atom2bin_multi;
    delete [] binhead_multi;
    delete [] bins_multi;

    // Realloacte at updated maxtypes

    maxgroups = n_multi_groups;

    atom2bin_multi = new int*[maxgroups]();
    binhead_multi = new int*[maxgroups]();
    bins_multi = new int*[maxgroups]();

    memory->destroy(nbinx_multi);
    memory->destroy(nbiny_multi);
    memory->destroy(nbinz_multi);
    memory->create(nbinx_multi, maxgroups, "neigh:nbinx_multi");
    memory->create(nbiny_multi, maxgroups, "neigh:nbiny_multi");
    memory->create(nbinz_multi, maxgroups, "neigh:nbinz_multi");

    memory->destroy(mbins_multi);
    memory->destroy(mbinx_multi);
    memory->destroy(mbiny_multi);
    memory->destroy(mbinz_multi);
    memory->create(mbins_multi, maxgroups, "neigh:mbins_multi");
    memory->create(mbinx_multi, maxgroups, "neigh:mbinx_multi");
    memory->create(mbiny_multi, maxgroups, "neigh:mbiny_multi");
    memory->create(mbinz_multi, maxgroups, "neigh:mbinz_multi");

    memory->destroy(mbinxlo_multi);
    memory->destroy(mbinylo_multi);
    memory->destroy(mbinzlo_multi);
    memory->create(mbinxlo_multi, maxgroups, "neigh:mbinxlo_multi");
    memory->create(mbinylo_multi, maxgroups, "neigh:mbinylo_multi");
    memory->create(mbinzlo_multi, maxgroups, "neigh:mbinzlo_multi");

    memory->destroy(binsizex_multi);
    memory->destroy(binsizey_multi);
    memory->destroy(binsizez_multi);
    memory->create(binsizex_multi, maxgroups, "neigh:binsizex_multi");
    memory->create(binsizey_multi, maxgroups, "neigh:binsizey_multi");
    memory->create(binsizez_multi, maxgroups, "neigh:binsizez_multi");

    memory->destroy(bininvx_multi);
    memory->destroy(bininvy_multi);
    memory->destroy(bininvz_multi);
    memory->create(bininvx_multi, maxgroups, "neigh:bininvx_multi");
    memory->create(bininvy_multi, maxgroups, "neigh:bininvy_multi");
    memory->create(bininvz_multi, maxgroups, "neigh:bininvz_multi");

    memory->destroy(maxbins_multi);
    memory->create(maxbins_multi, maxgroups, "neigh:maxbins_multi");
    
    // make sure reallocation occurs in bin_atoms_setup()
    for (n = 0; n < maxgroups; n++) {
      maxbins_multi[n] = 0;
    }
    maxatom = 0;
  }

  igroupmin = igroup_min();

  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3],bsubboxlo[3],bsubboxhi[3];
  double *cutghost = comm->cutghost;

  if (triclinic == 0) {
    bsubboxlo[0] = domain->sublo[0] - cutghost[0];
    bsubboxlo[1] = domain->sublo[1] - cutghost[1];
    bsubboxlo[2] = domain->sublo[2] - cutghost[2];
    bsubboxhi[0] = domain->subhi[0] + cutghost[0];
    bsubboxhi[1] = domain->subhi[1] + cutghost[1];
    bsubboxhi[2] = domain->subhi[2] + cutghost[2];
  } else {
    double lo[3],hi[3];
    lo[0] = domain->sublo_lamda[0] - cutghost[0];
    lo[1] = domain->sublo_lamda[1] - cutghost[1];
    lo[2] = domain->sublo_lamda[2] - cutghost[2];
    hi[0] = domain->subhi_lamda[0] + cutghost[0];
    hi[1] = domain->subhi_lamda[1] + cutghost[1];
    hi[2] = domain->subhi_lamda[2] + cutghost[2];
    domain->bbox(lo,hi,bsubboxlo,bsubboxhi);
  }

  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // For each grouping...

  for (n = 0; n < maxgroups; n++) {
    // binsize_user only relates to smallest type
    // optimal bin size is roughly 1/2 the type-type cutoff
    // special case of all cutoffs = 0.0, binsize = box size

    double binsize_optimal;
    if (n == igroupmin && binsizeflag) binsize_optimal = binsize_user;
    else binsize_optimal = 0.5*sqrt(cutmultisq[n][n]);
    if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
    double binsizeinv = 1.0/binsize_optimal;

    // test for too many global bins in any dimension due to huge global domain

    if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
	    bbox[2]*binsizeinv > MAXSMALLINT)
      error->all(FLERR,"Domain too large for neighbor bins");

    // create actual bins
    // always have one bin even if cutoff > bbox
    // for 2d, nbinz_multi[n] = 1

    nbinx_multi[n] = static_cast<int> (bbox[0]*binsizeinv);
    nbiny_multi[n] = static_cast<int> (bbox[1]*binsizeinv);
    if (dimension == 3) nbinz_multi[n] = static_cast<int> (bbox[2]*binsizeinv);
    else nbinz_multi[n] = 1;

    if (nbinx_multi[n] == 0) nbinx_multi[n] = 1;
    if (nbiny_multi[n] == 0) nbiny_multi[n] = 1;
    if (nbinz_multi[n] == 0) nbinz_multi[n] = 1;

    // compute actual bin size for nbins to fit into box exactly
    // error if actual bin size << cutoff, since will create a zillion bins
    // this happens when nbin = 1 and box size << cutoff
    // typically due to non-periodic, flat system in a particular dim
    // in that extreme case, should use NSQ not BIN neighbor style

    binsizex_multi[n] = bbox[0]/nbinx_multi[n];
    binsizey_multi[n] = bbox[1]/nbiny_multi[n];
    binsizez_multi[n] = bbox[2]/nbinz_multi[n];

    bininvx_multi[n] = 1.0 / binsizex_multi[n];
    bininvy_multi[n] = 1.0 / binsizey_multi[n];
    bininvz_multi[n] = 1.0 / binsizez_multi[n];

    if (binsize_optimal*bininvx_multi[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvy_multi[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvz_multi[n] > CUT2BIN_RATIO)
      error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

    // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
    // coord = lowest and highest values of coords for my ghost atoms
    // static_cast(-1.5) = -1, so subract additional -1
    // add in SMALL for round-off safety

    int mbinxhi,mbinyhi,mbinzhi;
    double coord;

    coord = bsubboxlo[0] - SMALL*bbox[0];
    mbinxlo_multi[n] = static_cast<int> ((coord-bboxlo[0])*bininvx_multi[n]);
    if (coord < bboxlo[0]) mbinxlo_multi[n] = mbinxlo_multi[n] - 1;
    coord = bsubboxhi[0] + SMALL*bbox[0];
    mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx_multi[n]);

    coord = bsubboxlo[1] - SMALL*bbox[1];
    mbinylo_multi[n] = static_cast<int> ((coord-bboxlo[1])*bininvy_multi[n]);
    if (coord < bboxlo[1]) mbinylo_multi[n] = mbinylo_multi[n] - 1;
    coord = bsubboxhi[1] + SMALL*bbox[1];
    mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy_multi[n]);

    if (dimension == 3) {
      coord = bsubboxlo[2] - SMALL*bbox[2];
      mbinzlo_multi[n] = static_cast<int> ((coord-bboxlo[2])*bininvz_multi[n]);
      if (coord < bboxlo[2]) mbinzlo_multi[n] = mbinzlo_multi[n] - 1;
      coord = bsubboxhi[2] + SMALL*bbox[2];
      mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz_multi[n]);
    }

    // extend bins by 1 to insure stencil extent is included
    // for 2d, only 1 bin in z

    mbinxlo_multi[n] = mbinxlo_multi[n] - 1;
    mbinxhi = mbinxhi + 1;
    mbinx_multi[n] = mbinxhi - mbinxlo_multi[n] + 1;

    mbinylo_multi[n] = mbinylo_multi[n] - 1;
    mbinyhi = mbinyhi + 1;
    mbiny_multi[n] = mbinyhi - mbinylo_multi[n] + 1;

    if (dimension == 3) {
      mbinzlo_multi[n] = mbinzlo_multi[n] - 1;
      mbinzhi = mbinzhi + 1;
    } else mbinzlo_multi[n] = mbinzhi = 0;
    mbinz_multi[n] = mbinzhi - mbinzlo_multi[n] + 1;

    bigint bbin = ((bigint) mbinx_multi[n]) 
      * ((bigint) mbiny_multi[n]) * ((bigint) mbinz_multi[n]) + 1;
    if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
    mbins_multi[n] = bbin;
  }

}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms by type
------------------------------------------------------------------------- */

void NBinMulti::bin_atoms()
{
  int i,ibin,n;

  last_bin = update->ntimestep;
  for (n = 0; n < maxgroups; n++) {
    for (i = 0; i < mbins_multi[n]; i++) binhead_multi[n][i] = -1;
  }

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        n = map_type_multi[type[i]];
        ibin = coord2bin_multi(x[i], n);
        atom2bin_multi[n][i] = ibin;
        bins_multi[n][i] = binhead_multi[n][ibin];
        binhead_multi[n][ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      n = map_type_multi[type[i]];
      ibin = coord2bin_multi(x[i], n);
      atom2bin_multi[n][i] = ibin;
      bins_multi[n][i] = binhead_multi[n][ibin];
      binhead_multi[n][ibin] = i;
    }
  } else {
    for (i = nall-1; i >= 0; i--) {
      n = map_type_multi[type[i]];
      ibin = coord2bin_multi(x[i], n);
      atom2bin_multi[n][i] = ibin;
      bins_multi[n][i] = binhead_multi[n][ibin];
      binhead_multi[n][ibin] = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

double NBinMulti::memory_usage()
{
  double bytes = 0;
  
  for (int m = 0; m < maxgroups; m++) {
    bytes += maxbins_multi[m]*sizeof(int);
    bytes += 2*maxatom*sizeof(int);
  }
  return bytes;
}
