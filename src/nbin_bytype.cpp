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

#include "nbin_bytype.h"
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

NBinBytype::NBinBytype(LAMMPS *lmp) : NBin(lmp) {

  nbinx_type = NULL; nbiny_type = NULL; nbinz_type = NULL;
  mbins_type = NULL;
  mbinx_type = NULL; mbiny_type = NULL, mbinz_type = NULL;
  mbinxlo_type = NULL;
  mbinylo_type = NULL;
  mbinzlo_type = NULL;
  binsizex_type = NULL; binsizey_type = NULL; binsizez_type = NULL;
  bininvx_type = NULL; bininvy_type = NULL; bininvz_type = NULL;
  binhead_type = NULL;
  bins_type = NULL;
  atom2bin_type = NULL;
  maxtypes = 0;
  maxbins_type = NULL;
}

NBinBytype::~NBinBytype() {

  memory->destroy(nbinx_type);
  memory->destroy(nbiny_type);
  memory->destroy(nbinz_type);
  memory->destroy(mbins_type);
  memory->destroy(mbinx_type);
  memory->destroy(mbiny_type);
  memory->destroy(mbinz_type);
  memory->destroy(mbinxlo_type);
  memory->destroy(mbinylo_type);
  memory->destroy(mbinzlo_type);

  memory->destroy(binsizex_type);
  memory->destroy(binsizey_type);
  memory->destroy(binsizez_type);
  memory->destroy(bininvx_type);
  memory->destroy(bininvy_type);
  memory->destroy(bininvz_type);

  for (int n = 1; n <= maxtypes; n++) {
    memory->destroy(binhead_type[n]);
    memory->destroy(bins_type[n]);
    memory->destroy(atom2bin_type[n]);
  }
  delete [] binhead_type;
  delete [] bins_type;
  delete [] atom2bin_type;

  memory->destroy(maxbins_type);
}

/* ----------------------------------------------------------------------
   arrange storage for types
   allows ntypes to change, but not currently expected after initialisation
   ------------------------------------------------------------------------- */

void NBinBytype::setup_types() {

  int n;

  if (atom->ntypes > maxtypes) {

    // Clear any/all memory for existing types

    for (n = 1; n <= maxtypes; n++) {
      memory->destroy(atom2bin_type[n]);
      memory->destroy(binhead_type[n]);
      memory->destroy(bins_type[n]);
    }
    delete [] atom2bin_type;
    delete [] binhead_type;
    delete [] bins_type;

    // Realloacte at updated maxtypes

    maxtypes = atom->ntypes;

    atom2bin_type = new int*[maxtypes+1]();
    binhead_type = new int*[maxtypes+1]();
    bins_type = new int*[maxtypes+1]();

    memory->destroy(nbinx_type);
    memory->destroy(nbiny_type);
    memory->destroy(nbinz_type);
    memory->create(nbinx_type, maxtypes+1, "nBinType:nbinx_type");
    memory->create(nbiny_type, maxtypes+1, "nBinType:nbiny_type");
    memory->create(nbinz_type, maxtypes+1, "nBinType:nbinz_type");

    memory->destroy(mbins_type);
    memory->destroy(mbinx_type);
    memory->destroy(mbiny_type);
    memory->destroy(mbinz_type);
    memory->create(mbins_type, maxtypes+1, "nBinType:mbins_type");
    memory->create(mbinx_type, maxtypes+1, "nBinType:mbinx_type");
    memory->create(mbiny_type, maxtypes+1, "nBinType:mbiny_type");
    memory->create(mbinz_type, maxtypes+1, "nBinType:mbinz_type");

    memory->destroy(mbinxlo_type);
    memory->destroy(mbinylo_type);
    memory->destroy(mbinzlo_type);
    memory->create(mbinxlo_type, maxtypes+1, "nBinType:mbinxlo_type");
    memory->create(mbinylo_type, maxtypes+1, "nBinType:mbinylo_type");
    memory->create(mbinzlo_type, maxtypes+1, "nBinType:mbinzlo_type");

    memory->destroy(binsizex_type);
    memory->destroy(binsizey_type);
    memory->destroy(binsizez_type);
    memory->create(binsizex_type, maxtypes+1, "nBinType:binsizex_type");
    memory->create(binsizey_type, maxtypes+1, "nBinType:binsizey_type");
    memory->create(binsizez_type, maxtypes+1, "nBinType:binsizez_type");

    memory->destroy(bininvx_type);
    memory->destroy(bininvy_type);
    memory->destroy(bininvz_type);
    memory->create(bininvx_type, maxtypes+1, "nBinType:bininvx_type");
    memory->create(bininvy_type, maxtypes+1, "nBinType:bininvy_type");
    memory->create(bininvz_type, maxtypes+1, "nBinType:bininvz_type");

    memory->destroy(maxbins_type);
    memory->create(maxbins_type, maxtypes+1, "nBinType:maxbins_type");
    // make sure reallocation occurs in bin_atoms_setup()
    for (n = 1; n <= maxtypes; n++) {
      maxbins_type[n] = 0;
    }
    maxatom = 0;
  }

}

/* ---------------------------------------------------------------------
   identify index of type with smallest cutoff
   --------------------------------------------------------------------- */

int NBinBytype::itype_min() {

  int itypemin = 1;
  double ** cutneighsq;

  cutneighsq = neighbor->cutneighsq;

  for (int n = 1; n <= atom->ntypes; n++) {
    if (cutneighsq[n][n] < cutneighsq[itypemin][itypemin]) {
      itypemin = n;
    }
  }

  return itypemin;
}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   ---------------------------------------------------------------------- */

void NBinBytype::setup_bins(int style)
{
  int n;
  int itypemin;

  setup_types();
  itypemin = itype_min();

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

  // For each type...

  for (n = 1; n <= atom->ntypes; n++) {
    // binsize_user only relates to smallest type
    // optimal bin size is roughly 1/2 the type-type cutoff
    // special case of all cutoffs = 0.0, binsize = box size

    double binsize_optimal;
    if (n == itypemin && binsizeflag) binsize_optimal = binsize_user;
    else binsize_optimal = 0.5*sqrt(neighbor->cutneighsq[n][n]);
    if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
    double binsizeinv = 1.0/binsize_optimal;

    // test for too many global bins in any dimension due to huge global domain

    if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
	    bbox[2]*binsizeinv > MAXSMALLINT)
      error->all(FLERR,"Domain too large for neighbor bins");

    // create actual bins
    // always have one bin even if cutoff > bbox
    // for 2d, nbinz_type[n] = 1

    nbinx_type[n] = static_cast<int> (bbox[0]*binsizeinv);
    nbiny_type[n] = static_cast<int> (bbox[1]*binsizeinv);
    if (dimension == 3) nbinz_type[n] = static_cast<int> (bbox[2]*binsizeinv);
    else nbinz_type[n] = 1;

    if (nbinx_type[n] == 0) nbinx_type[n] = 1;
    if (nbiny_type[n] == 0) nbiny_type[n] = 1;
    if (nbinz_type[n] == 0) nbinz_type[n] = 1;

    // compute actual bin size for nbins to fit into box exactly
    // error if actual bin size << cutoff, since will create a zillion bins
    // this happens when nbin = 1 and box size << cutoff
    // typically due to non-periodic, flat system in a particular dim
    // in that extreme case, should use NSQ not BIN neighbor style

    binsizex_type[n] = bbox[0]/nbinx_type[n];
    binsizey_type[n] = bbox[1]/nbiny_type[n];
    binsizez_type[n] = bbox[2]/nbinz_type[n];

    bininvx_type[n] = 1.0 / binsizex_type[n];
    bininvy_type[n] = 1.0 / binsizey_type[n];
    bininvz_type[n] = 1.0 / binsizez_type[n];

    if (binsize_optimal*bininvx_type[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvy_type[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvz_type[n] > CUT2BIN_RATIO)
      error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

    // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
    // coord = lowest and highest values of coords for my ghost atoms
    // static_cast(-1.5) = -1, so subract additional -1
    // add in SMALL for round-off safety

    int mbinxhi,mbinyhi,mbinzhi;
    double coord;

    coord = bsubboxlo[0] - SMALL*bbox[0];
    mbinxlo_type[n] = static_cast<int> ((coord-bboxlo[0])*bininvx_type[n]);
    if (coord < bboxlo[0]) mbinxlo_type[n] = mbinxlo_type[n] - 1;
    coord = bsubboxhi[0] + SMALL*bbox[0];
    mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx_type[n]);

    coord = bsubboxlo[1] - SMALL*bbox[1];
    mbinylo_type[n] = static_cast<int> ((coord-bboxlo[1])*bininvy_type[n]);
    if (coord < bboxlo[1]) mbinylo_type[n] = mbinylo_type[n] - 1;
    coord = bsubboxhi[1] + SMALL*bbox[1];
    mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy_type[n]);

    if (dimension == 3) {
      coord = bsubboxlo[2] - SMALL*bbox[2];
      mbinzlo_type[n] = static_cast<int> ((coord-bboxlo[2])*bininvz_type[n]);
      if (coord < bboxlo[2]) mbinzlo_type[n] = mbinzlo_type[n] - 1;
      coord = bsubboxhi[2] + SMALL*bbox[2];
      mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz_type[n]);
    }

    // extend bins by 1 to insure stencil extent is included
    // for 2d, only 1 bin in z

    mbinxlo_type[n] = mbinxlo_type[n] - 1;
    mbinxhi = mbinxhi + 1;
    mbinx_type[n] = mbinxhi - mbinxlo_type[n] + 1;

    mbinylo_type[n] = mbinylo_type[n] - 1;
    mbinyhi = mbinyhi + 1;
    mbiny_type[n] = mbinyhi - mbinylo_type[n] + 1;

    if (dimension == 3) {
      mbinzlo_type[n] = mbinzlo_type[n] - 1;
      mbinzhi = mbinzhi + 1;
    } else mbinzlo_type[n] = mbinzhi = 0;
    mbinz_type[n] = mbinzhi - mbinzlo_type[n] + 1;

    bigint bbin = ((bigint) mbinx_type[n]) 
      * ((bigint) mbiny_type[n]) * ((bigint) mbinz_type[n]) + 1;
    if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
    mbins_type[n] = bbin;
  }

}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms by type
------------------------------------------------------------------------- */

void NBinBytype::bin_atoms()
{
  int i,ibin,n;

  last_bin = update->ntimestep;
  for (n = 1; n <= maxtypes; n++) {
    for (i = 0; i < mbins_type[n]; i++) binhead_type[n][i] = -1;
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
        n = type[i];
        ibin = coord2bin(x[i], n);
        atom2bin_type[n][i] = ibin;
        bins_type[n][i] = binhead_type[n][ibin];
        binhead_type[n][ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      n = type[i];
      ibin = coord2bin(x[i], n);
      atom2bin_type[n][i] = ibin;
      bins_type[n][i] = binhead_type[n][ibin];
      binhead_type[n][ibin] = i;
    }
  } else {
    for (i = nall-1; i >= 0; i--) {
      n = type[i];
      ibin = coord2bin(x[i], n);
      atom2bin_type[n][i] = ibin;
      bins_type[n][i] = binhead_type[n][ibin];
      binhead_type[n][ibin] = i;
    }
  }
}


/* ----------------------------------------------------------------------
   allocate/reallocate vectors
------------------------------------------------------------------------- */

void NBinBytype::bin_atoms_setup(int nall) {

  // all atom2bin, bins must be of length nall
  if (nall > maxatom) {
    for (int n = 1; n <= maxtypes; n++) {
      memory->destroy(bins_type[n]);
      memory->destroy(atom2bin_type[n]);
      memory->create(bins_type[n], nall, "NBinBytype:bin_type");
      memory->create(atom2bin_type[n], nall, "NBinBytype:atom2bin_type");
    }
    maxatom = nall;
  }

  for (int n = 1; n <= maxtypes; n++) {
    if (mbins_type[n] > maxbins_type[n]) {
      maxbins_type[n] = mbins_type[n];
      memory->destroy(binhead_type[n]);
      memory->create(binhead_type[n], mbins_type[n], "NBinBytype:mbins_type");
    }
  }

}


/* ----------------------------------------------------------------------
   convert atom coords into local bin # for bin type it
------------------------------------------------------------------------- */

int NBinBytype::coord2bin(double *x, int it)
{
  int ix,iy,iz;
  int ibin;

  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2]))
    error->one(FLERR,"Non-numeric positions - simulation unstable");

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx_type[it]) + nbinx_type[it];
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_type[it]);
    ix = MIN(ix,nbinx_type[it]-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx_type[it]) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy_type[it]) + nbiny_type[it];
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_type[it]);
    iy = MIN(iy,nbiny_type[it]-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy_type[it]) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz_type[it]) + nbinz_type[it];
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_type[it]);
    iz = MIN(iz,nbinz_type[it]-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz_type[it]) - 1;

  
  ibin = (iz-mbinzlo_type[it])*mbiny_type[it]*mbinx_type[it]
       + (iy-mbinylo_type[it])*mbinx_type[it]
       + (ix-mbinxlo_type[it]);
  return ibin;
}


/* ----------------------------------------------------------------------
   tot up for types
   ---------------------------------------------------------------------- */

bigint NBinBytype::memory_usage()
{
  bigint bytes = 0;

  for (int m = 1; m < maxtypes; m++) {
    bytes += 2*maxatom*sizeof(int);
    bytes += maxbins_type[m]*sizeof(int);
  }
  return bytes;
}
