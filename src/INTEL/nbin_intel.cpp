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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "nbin_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NBinIntel::NBinIntel(LAMMPS *lmp) : NBinStandard(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
  _precision_mode = _fix->precision();
  _atombin = nullptr;
  _binpacked = nullptr;
}

/* ---------------------------------------------------------------------- */

NBinIntel::~NBinIntel() {
  memory->destroy(_atombin);
  memory->destroy(_binpacked);
}

/* ----------------------------------------------------------------------
   setup for bin_atoms()
------------------------------------------------------------------------- */

void NBinIntel::bin_atoms_setup(int nall)
{
  // binhead = per-bin vector, mbins in length
  // add 1 bin for INTEL package

  if (mbins > maxbin) {

    maxbin = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxbin+1,"neigh:binhead");

  }

  // bins = per-atom vector

  if (nall > maxatom) {
    maxatom = nall;

    memory->destroy(bins);
    memory->destroy(_atombin);
    memory->destroy(_binpacked);

    memory->create(bins,maxatom,"neigh:bins");
    memory->create(_atombin,maxatom,"neigh:bins");
    memory->create(_binpacked,maxatom,"neigh:bins");

    if (_precision_mode == FixIntel::PREC_MODE_MIXED)
      _fix->get_mixed_buffers()->set_bininfo(_atombin,_binpacked);
    else if (_precision_mode == FixIntel::PREC_MODE_SINGLE)
      _fix->get_single_buffers()->set_bininfo(_atombin,_binpacked);
    else
      _fix->get_double_buffers()->set_bininfo(_atombin,_binpacked);
  }
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

  // ---------- Sanity check for padding --------------
  {
    const flt_t dx = (INTEL_BIGP - bboxhi[0]);
    const flt_t dy = (INTEL_BIGP - bboxhi[1]);
    const flt_t dz = (INTEL_BIGP - bboxhi[2]);
    if (dx * dx + dy * dy + dz * dz < static_cast<flt_t>(neighbor->cutneighmaxsq))
      error->one(FLERR,"INTEL package expects no atoms within cutoff of (1e15,1e15,1e15).");
  }

  // ---------- Grow and cast/pack buffers -------------
  _fix->start_watch(TIME_PACK);
  buffers->grow(nall, atom->nlocal, comm->nthreads);

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
