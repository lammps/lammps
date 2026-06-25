// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "nbin_manual.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include <unordered_set>

using namespace LAMMPS_NS;

static constexpr double SMALL = 1.0e-6;
static constexpr double CUT2BIN_RATIO = 100.0;

/* ---------------------------------------------------------------------- */

NBinManual::NBinManual(LAMMPS *lmp) : NBin(lmp)
{
  cutoff_custom = 0.0;
  maxbin_custom = maxcustom = 0;
  binhead_custom = nullptr;
  bins_custom = nullptr;
  custom2bin = nullptr;
}

/* ----------------------------------------------------------------------
   manually assign neighbor info
     then calculate setup information
------------------------------------------------------------------------- */

void NBinManual::assign_neighbor_info(double cutmax)
{
  cutneighmax = cutmax;

  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;
}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinManual::bin_atoms_setup(int nall)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for INTEL package
  if (mbins > maxbin) {
    maxbin = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxbin,"neigh:binhead");
  }

  // bins and atom2bin = per-atom vectors
  // for both local and ghost atoms

  if (nall > maxatom) {
    maxatom = nall;
    memory->destroy(bins);
    memory->create(bins,maxatom,"neigh:bins");
    memory->destroy(atom2bin);
    memory->create(atom2bin,maxatom,"neigh:atom2bin");
  }
}

/* ----------------------------------------------------------------------
   analog of setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinManual::bin_custom_setup(double **xcustom, int ncustom)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for INTEL package

  if (mbins > maxbin_custom) {
    maxbin_custom = mbins;
    memory->destroy(binhead_custom);
    memory->create(binhead_custom,maxbin_custom,"neigh:binhead_custom");
  }

  int nall = ncustom;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int dim = domain->dimension;

  double limitlo[3], limithi[3];
  for (int a = 0; a < dim; a++) {
    limitlo[a] = boxlo[a] + cutneighmax;
    limithi[a] = boxhi[a] - cutneighmax;
  }

  for (int i = 0; i < ncustom; i++) {
    for (int a = 0; a < dim; a++) {
      if (!domain->periodicity[a]) continue;
      if (xcustom[i][a] < limitlo[a]) nall += 1;
      if (xcustom[i][a] > limithi[a]) nall += 1;
    }
  }

  // bins and atom2bin = per-atom vectors
  // for both local and ghost atoms

  if (nall > maxcustom) {
    maxcustom = nall;
    memory->destroy(bins_custom);
    memory->create(bins_custom,maxcustom,"neigh:bins_custom");
    memory->destroy(custom2bin);
    memory->create(custom2bin,maxcustom,"neigh:custom2bin");
  }
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

void NBinManual::setup_bins(int style)
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by cutghost or cutneighmax
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3],bsubboxlo[3],bsubboxhi[3];
  double *cutghost = comm->cutghost;

  // cutneighmax can be larger than cut ghost, so take max to calculate bins

  double cut_limit[3];
  if (triclinic == 0) {
    cut_limit[0] = MAX(cutneighmax, cutghost[0]);
    cut_limit[1] = MAX(cutneighmax, cutghost[1]);
    cut_limit[2] = MAX(cutneighmax, cutghost[2]);

    bsubboxlo[0] = domain->sublo[0] - cut_limit[0];
    bsubboxlo[1] = domain->sublo[1] - cut_limit[1];
    bsubboxlo[2] = domain->sublo[2] - cut_limit[2];
    bsubboxhi[0] = domain->subhi[0] + cut_limit[0];
    bsubboxhi[1] = domain->subhi[1] + cut_limit[1];
    bsubboxhi[2] = domain->subhi[2] + cut_limit[2];
  } else {
    double *h_inv = domain->h_inv;
    double length0,length1,length2;
    length0 = sqrt(h_inv[0]*h_inv[0] + h_inv[5]*h_inv[5] + h_inv[4]*h_inv[4]);
    length1 = sqrt(h_inv[1]*h_inv[1] + h_inv[3]*h_inv[3]);
    length2 = h_inv[2];

    cut_limit[0] = MAX(cutneighmax * length0, cutghost[0]);
    cut_limit[1] = MAX(cutneighmax * length1, cutghost[1]);
    cut_limit[2] = MAX(cutneighmax * length2, cutghost[2]);

    double lo[3],hi[3];
    lo[0] = domain->sublo_lamda[0] - cut_limit[0];
    lo[1] = domain->sublo_lamda[1] - cut_limit[1];
    lo[2] = domain->sublo_lamda[2] - cut_limit[2];
    hi[0] = domain->subhi_lamda[0] + cut_limit[0];
    hi[1] = domain->subhi_lamda[1] + cut_limit[1];
    hi[2] = domain->subhi_lamda[2] + cut_limit[2];
    domain->bbox(lo,hi,bsubboxlo,bsubboxhi);
  }

  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  double binsize_optimal;
  binsize_optimal = 0.5*cutneighmax;
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

  // extend bins by 1 to ensure stencil extent is included
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
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void NBinManual::bin_custom(double **xcustom, int ncustom)
{
  int i, ibin;

  last_bin = update->ntimestep;
  for (i = 0; i < mbins; i++) binhead_custom[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost objects (assumed to be at end of xcustom) at end of list

  for (i = ncustom-1; i >= 0; i--) {
    ibin = coord2bin(xcustom[i]);

    // check within bin limits (does not require for ghosts)
    if (ibin >= 0 && ibin < maxbin_custom) {
      custom2bin[i] = ibin;
      bins_custom[i] = binhead_custom[ibin];
      binhead_custom[ibin] = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void NBinManual::bin_atoms()
{
  int i,ibin;

  last_bin = update->ntimestep;
  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal; // + atom->nghost; // JTC I don't think I need these

  for (i = nall-1; i >= 0; i--) {
    ibin = coord2bin(x[i]);
    atom2bin[i] = ibin;
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }
}

/* ---------------------------------------------------------------------- */

double NBinManual::memory_usage()
{
  double bytes = 0;
  bytes += (double)maxbin*sizeof(int);
  bytes += (double)2*maxatom*sizeof(int);
  return bytes;
}
