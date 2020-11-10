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

#include "nbin_multi2.h"
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

NBinMulti2::NBinMulti2(LAMMPS *lmp) : NBin(lmp) {}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinMulti2::bin_atoms_setup(int nall) 
{
  // binhead_multi2[n] = per-bin vector mbins in length mbins_multi2[n]
  
  for (int n = 1; n <= maxtypes; n++) {
    if (mbins_multi2[n] > maxbins_multi2[n]) {
      maxbins_multi2[n] = mbins_multi2[n];
      memory->destroy(binhead_multi2[n]);
      memory->create(binhead_multi2[n], mbins_multi2[n], "neigh:mbins_multi2");
    }
  }

  // bins_multi2[n] and atom2bin_multi2[n] = per-atom vectors
  // for both local and ghost atoms

  if (nall > maxatom) {
    maxatom = nall;
    for (int n = 1; n <= maxtypes; n++) {
      memory->destroy(bins_multi2[n]);
      memory->destroy(atom2bin_multi2[n]);
      memory->create(bins_multi2[n], maxatom, "neigh:bin_multi2");
      memory->create(atom2bin_multi2[n], maxatom, "neigh:atom2bin_multi2");
    }
  }
}

/* ---------------------------------------------------------------------
   Identify index of type with smallest cutoff
------------------------------------------------------------------------ */

int NBinMulti2::itype_min() {

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

void NBinMulti2::setup_bins(int style)
{
  int n;
  int itypemin;

  // Initialize arrays
  
  if (atom->ntypes > maxtypes) {

    // Clear any/all memory for existing types

    for (n = 1; n <= maxtypes; n++) {
      memory->destroy(atom2bin_multi2[n]);
      memory->destroy(binhead_multi2[n]);
      memory->destroy(bins_multi2[n]);
    }
    delete [] atom2bin_multi2;
    delete [] binhead_multi2;
    delete [] bins_multi2;

    // Realloacte at updated maxtypes

    maxtypes = atom->ntypes;

    atom2bin_multi2 = new int*[maxtypes+1]();
    binhead_multi2 = new int*[maxtypes+1]();
    bins_multi2 = new int*[maxtypes+1]();

    memory->destroy(nbinx_multi2);
    memory->destroy(nbiny_multi2);
    memory->destroy(nbinz_multi2);
    memory->create(nbinx_multi2, maxtypes+1, "neigh:nbinx_multi2");
    memory->create(nbiny_multi2, maxtypes+1, "neigh:nbiny_multi2");
    memory->create(nbinz_multi2, maxtypes+1, "neigh:nbinz_multi2");

    memory->destroy(mbins_multi2);
    memory->destroy(mbinx_multi2);
    memory->destroy(mbiny_multi2);
    memory->destroy(mbinz_multi2);
    memory->create(mbins_multi2, maxtypes+1, "neigh:mbins_multi2");
    memory->create(mbinx_multi2, maxtypes+1, "neigh:mbinx_multi2");
    memory->create(mbiny_multi2, maxtypes+1, "neigh:mbiny_multi2");
    memory->create(mbinz_multi2, maxtypes+1, "neigh:mbinz_multi2");

    memory->destroy(mbinxlo_multi2);
    memory->destroy(mbinylo_multi2);
    memory->destroy(mbinzlo_multi2);
    memory->create(mbinxlo_multi2, maxtypes+1, "neigh:mbinxlo_multi2");
    memory->create(mbinylo_multi2, maxtypes+1, "neigh:mbinylo_multi2");
    memory->create(mbinzlo_multi2, maxtypes+1, "neigh:mbinzlo_multi2");

    memory->destroy(binsizex_multi2);
    memory->destroy(binsizey_multi2);
    memory->destroy(binsizez_multi2);
    memory->create(binsizex_multi2, maxtypes+1, "neigh:binsizex_multi2");
    memory->create(binsizey_multi2, maxtypes+1, "neigh:binsizey_multi2");
    memory->create(binsizez_multi2, maxtypes+1, "neigh:binsizez_multi2");

    memory->destroy(bininvx_multi2);
    memory->destroy(bininvy_multi2);
    memory->destroy(bininvz_multi2);
    memory->create(bininvx_multi2, maxtypes+1, "neigh:bininvx_multi2");
    memory->create(bininvy_multi2, maxtypes+1, "neigh:bininvy_multi2");
    memory->create(bininvz_multi2, maxtypes+1, "neigh:bininvz_multi2");

    memory->destroy(maxbins_multi2);
    memory->create(maxbins_multi2, maxtypes+1, "neigh:maxbins_multi2");
    // make sure reallocation occurs in bin_atoms_setup()
    for (n = 1; n <= maxtypes; n++) {
      maxbins_multi2[n] = 0;
    }
    maxatom = 0;
  }

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
    // for 2d, nbinz_multi2[n] = 1

    nbinx_multi2[n] = static_cast<int> (bbox[0]*binsizeinv);
    nbiny_multi2[n] = static_cast<int> (bbox[1]*binsizeinv);
    if (dimension == 3) nbinz_multi2[n] = static_cast<int> (bbox[2]*binsizeinv);
    else nbinz_multi2[n] = 1;

    if (nbinx_multi2[n] == 0) nbinx_multi2[n] = 1;
    if (nbiny_multi2[n] == 0) nbiny_multi2[n] = 1;
    if (nbinz_multi2[n] == 0) nbinz_multi2[n] = 1;

    // compute actual bin size for nbins to fit into box exactly
    // error if actual bin size << cutoff, since will create a zillion bins
    // this happens when nbin = 1 and box size << cutoff
    // typically due to non-periodic, flat system in a particular dim
    // in that extreme case, should use NSQ not BIN neighbor style

    binsizex_multi2[n] = bbox[0]/nbinx_multi2[n];
    binsizey_multi2[n] = bbox[1]/nbiny_multi2[n];
    binsizez_multi2[n] = bbox[2]/nbinz_multi2[n];

    bininvx_multi2[n] = 1.0 / binsizex_multi2[n];
    bininvy_multi2[n] = 1.0 / binsizey_multi2[n];
    bininvz_multi2[n] = 1.0 / binsizez_multi2[n];

    if (binsize_optimal*bininvx_multi2[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvy_multi2[n] > CUT2BIN_RATIO ||
	    binsize_optimal*bininvz_multi2[n] > CUT2BIN_RATIO)
      error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

    // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
    // coord = lowest and highest values of coords for my ghost atoms
    // static_cast(-1.5) = -1, so subract additional -1
    // add in SMALL for round-off safety

    int mbinxhi,mbinyhi,mbinzhi;
    double coord;

    coord = bsubboxlo[0] - SMALL*bbox[0];
    mbinxlo_multi2[n] = static_cast<int> ((coord-bboxlo[0])*bininvx_multi2[n]);
    if (coord < bboxlo[0]) mbinxlo_multi2[n] = mbinxlo_multi2[n] - 1;
    coord = bsubboxhi[0] + SMALL*bbox[0];
    mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx_multi2[n]);

    coord = bsubboxlo[1] - SMALL*bbox[1];
    mbinylo_multi2[n] = static_cast<int> ((coord-bboxlo[1])*bininvy_multi2[n]);
    if (coord < bboxlo[1]) mbinylo_multi2[n] = mbinylo_multi2[n] - 1;
    coord = bsubboxhi[1] + SMALL*bbox[1];
    mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy_multi2[n]);

    if (dimension == 3) {
      coord = bsubboxlo[2] - SMALL*bbox[2];
      mbinzlo_multi2[n] = static_cast<int> ((coord-bboxlo[2])*bininvz_multi2[n]);
      if (coord < bboxlo[2]) mbinzlo_multi2[n] = mbinzlo_multi2[n] - 1;
      coord = bsubboxhi[2] + SMALL*bbox[2];
      mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz_multi2[n]);
    }

    // extend bins by 1 to insure stencil extent is included
    // for 2d, only 1 bin in z

    mbinxlo_multi2[n] = mbinxlo_multi2[n] - 1;
    mbinxhi = mbinxhi + 1;
    mbinx_multi2[n] = mbinxhi - mbinxlo_multi2[n] + 1;

    mbinylo_multi2[n] = mbinylo_multi2[n] - 1;
    mbinyhi = mbinyhi + 1;
    mbiny_multi2[n] = mbinyhi - mbinylo_multi2[n] + 1;

    if (dimension == 3) {
      mbinzlo_multi2[n] = mbinzlo_multi2[n] - 1;
      mbinzhi = mbinzhi + 1;
    } else mbinzlo_multi2[n] = mbinzhi = 0;
    mbinz_multi2[n] = mbinzhi - mbinzlo_multi2[n] + 1;

    bigint bbin = ((bigint) mbinx_multi2[n]) 
      * ((bigint) mbiny_multi2[n]) * ((bigint) mbinz_multi2[n]) + 1;
    if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
    mbins_multi2[n] = bbin;
  }

}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms by type
------------------------------------------------------------------------- */

void NBinMulti2::bin_atoms()
{
  int i,ibin,n;

  last_bin = update->ntimestep;
  for (n = 1; n <= maxtypes; n++) {
    for (i = 0; i < mbins_multi2[n]; i++) binhead_multi2[n][i] = -1;
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
        atom2bin_multi2[n][i] = ibin;
        bins_multi2[n][i] = binhead_multi2[n][ibin];
        binhead_multi2[n][ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      n = type[i];
      ibin = coord2bin(x[i], n);
      atom2bin_multi2[n][i] = ibin;
      bins_multi2[n][i] = binhead_multi2[n][ibin];
      binhead_multi2[n][ibin] = i;
    }
  } else {
    for (i = nall-1; i >= 0; i--) {
      n = type[i];
      ibin = coord2bin(x[i], n);
      atom2bin_multi2[n][i] = ibin;
      bins_multi2[n][i] = binhead_multi2[n][ibin];
      binhead_multi2[n][ibin] = i;
    }
  }
}

/* ----------------------------------------------------------------------
   convert atom coords into local bin # for a particular type
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int NBinMulti2::coord2bin(double *x, int it)
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


/* ---------------------------------------------------------------------- */

double NBinMulti2::memory_usage()
{
  double bytes = 0;
  
  for (int m = 1; m < maxtypes; m++) {
    bytes += maxbins_multi2[m]*sizeof(int);
    bytes += 2*maxatom*sizeof(int);
  }
  return bytes;
}
