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

#include "pair_zero_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairZeroIntel::PairZeroIntel(LAMMPS *lmp) : PairZero(lmp), fix(nullptr)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairZeroIntel::compute(int eflag, int vflag)
{
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float, double>(eflag, vflag, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double, double>(eflag, vflag, fix->get_double_buffers());
  else
    compute<float, float>(eflag, vflag, fix->get_single_buffers());

  fix->balance_stamp();
  vflag_fdotr = 0;
}

/* ----------------------------------------------------------------------
   no pair interactions, but provide the packed coordinates and the
   zeroed force and energy buffers that the bond, angle, dihedral, and
   improper styles of the INTEL package add their contributions to
------------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairZeroIntel::compute(int eflag, int vflag, IntelBuffers<flt_t, acc_t> *buffers)
{
  ev_init(eflag, vflag);
  if (vflag_atom) error->all(FLERR, "INTEL package does not support per-atom stress");

  const int nthreads = comm->nthreads;
  const int ago = neighbor->ago;

  // on steps with a neighbor list build, the neighbor list build packs the coordinates

  if (ago != 0 && fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);

    int packthreads;
    if (nthreads > INTEL_HTHREADS)
      packthreads = nthreads;
    else
      packthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel if (packthreads > 1)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost, packthreads,
                                sizeof(ATOM_T));
      buffers->thr_pack(ifrom, ito, ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  int nlocal, nall, minlocal;
  fix->get_buffern(nlocal, nall, minlocal);

  int f_stride;
  IP_PRE_get_transfern(force->newton_pair, buffers, f_stride);

  int tc;
  FORCE_T *_noalias f_start;
  acc_t *_noalias ev_global;
  IP_PRE_get_buffers(buffers, fix, tc, f_start, ev_global);

  // with newton on, each thread has its own force buffer

  if (force->newton_pair)
    memset(f_start, 0, (size_t) tc * f_stride * sizeof(FORCE_T));
  else
    memset(f_start, 0, (size_t) f_stride * sizeof(FORCE_T));

  if (eflag || vflag)
    for (int i = 0; i < 8; ++i) ev_global[i] = (acc_t) 0;

  fix->stop_watch(TIME_HOST_PAIR);

  // no fdotr virial reduction was done here, so the per-thread force buffers must be reduced

  if (eflag || vflag)
    fix->add_result_array(f_start, ev_global, eflag_atom);
  else
    fix->add_result_array(f_start, nullptr);
}

/* ---------------------------------------------------------------------- */

void PairZeroIntel::init_style()
{
  PairZero::init_style();
  if (force->newton_pair == 0) neighbor->find_request(this)->enable_full();

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();
}
