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

#include "neigh_bin.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighBin::NeighBin(LAMMPS *lmp) : Pointers(lmp)
{
  last_setup = last_bin = last_bin_memory = -1;
  maxbin = maxatom = 0;
  binhead = NULL;
  bins = NULL;

  // geometry settings

  dimension = domain->dimension;
  triclinic = domain->triclinic;
}

/* ---------------------------------------------------------------------- */

NeighBin::~NeighBin()
{
  memory->destroy(binhead);
  memory->destroy(bins);
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class
------------------------------------------------------------------------- */

void NeighBin::copy_neighbor_info()
{
  includegroup = neighbor->includegroup;
  cutneighmin = neighbor->cutneighmin;
  cutneighmax = neighbor->cutneighmax;
  binsizeflag = neighbor->binsizeflag;
  binsize_user = neighbor->binsize_user;
  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;
}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */

enum{NSQ,BIN,MULTI};       // also in Neighbor

#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100

void NeighBin::setup_bins(int style)
{
  last_setup = update->ntimestep;

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

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else if (style == BIN) binsize_optimal = 0.5*cutneighmax;
  else binsize_optimal = 0.5*cutneighmin;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain

  if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
      bbox[2]*binsizeinv > MAXSMALLINT)
    error->all(FLERR,"Domain too large for neighbor bins");

  // create actual bins
  // always have one bin even if cutoff > bbox
  // for 2d, nbinz = 1

  nbinx = static_cast<int> (bbox[0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1]*binsizeinv);
  if (dimension == 3) nbinz = static_cast<int> (bbox[2]*binsizeinv);
  else nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly
  // error if actual bin size << cutoff, since will create a zillion bins
  // this happens when nbin = 1 and box size << cutoff
  // typically due to non-periodic, flat system in a particular dim
  // in that extreme case, should use NSQ not BIN neighbor style

  binsizex = bbox[0]/nbinx;
  binsizey = bbox[1]/nbiny;
  binsizez = bbox[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  if (binsize_optimal*bininvx > CUT2BIN_RATIO ||
      binsize_optimal*bininvy > CUT2BIN_RATIO ||
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  int mbinxhi,mbinyhi,mbinzhi;
  double coord;

  coord = bsubboxlo[0] - SMALL*bbox[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*bbox[0];
  mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*bbox[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*bbox[1];
  mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

  if (dimension == 3) {
    coord = bsubboxlo[2] - SMALL*bbox[2];
    mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
    if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
    coord = bsubboxhi[2] + SMALL*bbox[2];
    mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);
  }

  // extend bins by 1 to insure stencil extent is included
  // for 2d, only 1 bin in z

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  if (dimension == 3) {
    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;
  } else mbinzlo = mbinzhi = 0;
  mbinz = mbinzhi - mbinzlo + 1;

  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz) + 1;
  if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
  mbins = bbin;
}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NeighBin::bin_atoms_setup(int nall)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for USER-INTEL package

  if (mbins > maxbin) {
    maxbin = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxbin,"neigh:binhead");
    last_bin_memory = update->ntimestep;
  }

  // bins = per-atom vector

  if (nall > maxatom) {
    maxatom = nall;
    memory->destroy(bins);
    memory->create(bins,maxatom,"neigh:bins");
    last_bin_memory = update->ntimestep;
  }

  last_bin = update->ntimestep;
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

int NeighBin::coord2bin(double *x)
{
  int ix,iy,iz;

  if (!ISFINITE(x[0]) || !ISFINITE(x[1]) || !ISFINITE(x[2]))
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

/* ---------------------------------------------------------------------- */

bigint NeighBin::memory_usage()
{
  bigint bytes = 0;
  bytes += maxbin*sizeof(int);
  bytes += maxatom*sizeof(int);
  return bytes;
}
