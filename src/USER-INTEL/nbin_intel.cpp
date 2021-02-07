/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "nbin_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

#define SMALL 1.0e-6
#define CUT2BIN_RATIO 100

/* ---------------------------------------------------------------------- */

NBinIntel::NBinIntel(LAMMPS *lmp) : NBinStandard(lmp) {
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  _fix = static_cast<FixIntel *>(modify->fix[ifix]);
  _precision_mode = _fix->precision();
  _atombin = nullptr;
  _binpacked = nullptr;
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = _fix->coprocessor_number();
  _offload_alloc = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

NBinIntel::~NBinIntel() {
  #ifdef _LMP_INTEL_OFFLOAD
  if (_offload_alloc) {
    const int * binhead = this->binhead;
    const int * bins = this->bins;
    const int * _atombin = this->_atombin;
    const int * _binpacked = this->_binpacked;
    #pragma offload_transfer target(mic:_cop)   \
      nocopy(binhead,bins,_atombin,_binpacked:alloc_if(0) free_if(1))
  }
  #endif
  memory->destroy(_atombin);
  memory->destroy(_binpacked);
}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinIntel::bin_atoms_setup(int nall)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for USER-INTEL package

  if (mbins > maxbin) {
    #ifdef _LMP_INTEL_OFFLOAD
    if (_offload_alloc) {
      const int * binhead = this->binhead;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(binhead:alloc_if(0) free_if(1))
    }
    #endif

    maxbin = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxbin+1,"neigh:binhead");

    #ifdef _LMP_INTEL_OFFLOAD
    if (_fix->offload_balance() != 0) {
      int * binhead = this->binhead;
      #pragma offload_transfer target(mic:_cop) \
         nocopy(binhead:length(maxbin+1) alloc_if(1) free_if(0))
    }
    #endif
  }

  // bins = per-atom vector

  if (nall > maxatom) {
    maxatom = nall;

    #ifdef _LMP_INTEL_OFFLOAD
    if (_offload_alloc) {
      const int * bins = this->bins;
      const int * _atombin = this->_atombin;
      const int * _binpacked = this->_binpacked;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(bins,_atombin,_binpacked:alloc_if(0) free_if(1))
    }
    #endif
    memory->destroy(bins);
    memory->destroy(_atombin);
    memory->destroy(_binpacked);

    memory->create(bins,maxatom,"neigh:bins");
    memory->create(_atombin,maxatom,"neigh:bins");
    memory->create(_binpacked,maxatom,"neigh:bins");
    #ifdef _LMP_INTEL_OFFLOAD
    if (_fix->offload_balance() != 0) {
      const int * bins = this->bins;
      const int * _atombin = this->_atombin;
      const int * _binpacked = this->_binpacked;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(bins,_atombin,_binpacked:length(maxatom) alloc_if(1) free_if(0))
      _offload_alloc=1;
    }
    #endif

    if (_precision_mode == FixIntel::PREC_MODE_MIXED)
      _fix->get_mixed_buffers()->set_bininfo(_atombin,_binpacked);
    else if (_precision_mode == FixIntel::PREC_MODE_SINGLE)
      _fix->get_single_buffers()->set_bininfo(_atombin,_binpacked);
    else
      _fix->get_double_buffers()->set_bininfo(_atombin,_binpacked);
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

void NBinIntel::setup_bins(int style)
{
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
  // for MULTI_OLD style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else if (style == Neighbor::BIN) binsize_optimal = 0.5*cutneighmax;
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
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void NBinIntel::bin_atoms()
{
  last_bin = update->ntimestep;

  if (_precision_mode == FixIntel::PREC_MODE_MIXED)
    bin_atoms(_fix->get_mixed_buffers());
  else if (_precision_mode == FixIntel::PREC_MODE_SINGLE)
    bin_atoms(_fix->get_single_buffers());
  else
    bin_atoms(_fix->get_double_buffers());
}

template <class flt_t, class acc_t>
void NBinIntel::bin_atoms(IntelBuffers<flt_t,acc_t> * buffers) {
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int aend = _fix->offload_end_neighbor();


  // ---------- Sanity check for padding --------------
  {
    const flt_t dx = (INTEL_BIGP - bboxhi[0]);
    const flt_t dy = (INTEL_BIGP - bboxhi[1]);
    const flt_t dz = (INTEL_BIGP - bboxhi[2]);
    if (dx * dx + dy * dy + dz * dz <
        static_cast<flt_t>(neighbor->cutneighmaxsq))
      error->one(FLERR,
        "Intel package expects no atoms within cutoff of {1e15,1e15,1e15}.");
  }

  // ---------- Grow and cast/pack buffers -------------
  _fix->start_watch(TIME_PACK);
  buffers->grow(nall, atom->nlocal, comm->nthreads, aend);

  ATOM_T biga;
  biga.x = INTEL_BIGP;
  biga.y = INTEL_BIGP;
  biga.z = INTEL_BIGP;
  biga.w = 1;
  buffers->get_x()[nall] = biga;

  int nthreads;
  if (comm->nthreads > INTEL_HTHREADS) nthreads = comm->nthreads;
  else nthreads = 1;
  #if defined(_OPENMP)
  #pragma omp parallel if (nthreads > INTEL_HTHREADS)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id_align(ifrom, ito, tid, nall, nthreads,
                              sizeof(ATOM_T));
    buffers->thr_pack(ifrom, ito, 0);
  }
  _fix->stop_watch(TIME_PACK);


  // ---------- Bin Atoms -------------
  _fix->start_watch(TIME_HOST_NEIGHBOR);
  int * _noalias const atombin = this->_atombin;
  int * _noalias const binpacked = this->_binpacked;

  int i, ibin;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  int *mask = atom->mask;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        ibin = coord2bin(atom->x[i]);
        // Only necessary to store when neighboring ghost
        atombin[i] = ibin;
        bins[i] = binhead[ibin];
        binhead[ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      ibin = coord2bin(atom->x[i]);
      atombin[i] = ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  } else {
    for (i = nall-1; i >= 0; i--) {
      ibin = coord2bin(atom->x[i]);
      // Only necessary to store for ghost when neighboring ghost
      atombin[i] = ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  }
  int newhead = 0;
  for (i = 0; i < mbins; i++) {
    int j = binhead[i];
    binhead[i] = newhead;
    for ( ; j >= 0; j = bins[j])
      binpacked[newhead++] = j;
  }
  binhead[mbins] = newhead;
}

/* ---------------------------------------------------------------------- */

double NBinIntel::memory_usage()
{
  return NBinStandard::memory_usage() + maxatom*2*sizeof(int);
}
