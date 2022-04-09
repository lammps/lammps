/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "electrode_accel_intel.h"

#include "comm.h"
#include "fix_intel.h"
#include "intel_buffers.h"
#include "modify.h"

using namespace LAMMPS_NS;

ElectrodeAccelIntel::ElectrodeAccelIntel(class LAMMPS *lmp) : ElectrodeAccelInterface(lmp) {}

void ElectrodeAccelIntel::intel_find_fix()
{
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0) error->all(FLERR, "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);
}

void ElectrodeAccelIntel::intel_pack_buffers()
{
  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      intel_pack_buffers_prec<float, double>(fix->get_mixed_buffers());
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      intel_pack_buffers_prec<double, double>(fix->get_double_buffers());
      break;
    default:    // FixIntel::PREC_MODE_SINGLE
      intel_pack_buffers_prec<float, float>(fix->get_single_buffers());
  }
}

template <class flt_t, class acc_t>
void ElectrodeAccelIntel::intel_pack_buffers_prec(IntelBuffers<flt_t, acc_t> *buffers)
{
  fix->start_watch(TIME_PACK);
  int packthreads;
  if (comm->nthreads > INTEL_HTHREADS)
    packthreads = comm->nthreads;
  else
    packthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel if (packthreads > 1)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost, packthreads,
                              sizeof(ATOM_T));
    buffers->thr_pack(ifrom, ito, 0);
  }
  fix->stop_watch(TIME_PACK);
}
