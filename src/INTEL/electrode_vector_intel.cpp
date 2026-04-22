/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Shern Tee (UQ)
 ------------------------------------------------------------------------- */

#include "electrode_vector_intel.h"

#include "atom.h"
#include "comm.h"
#include "electrode_math.h"
#include "force.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cassert>

using namespace LAMMPS_NS;

void ElectrodeVectorIntel::get_fix_intel()
{
  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
}

void ElectrodeVectorIntel::pair_contribution(double *vector)
{
  _lrt = fix->lrt();
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
      pair_contribution<float,double>(fix->get_mixed_buffers(), vector);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      pair_contribution<double,double>(fix->get_double_buffers(), vector);
  else if (fix->precision() == FixIntel::PREC_MODE_SINGLE)
      pair_contribution<float,float>(fix->get_single_buffers(), vector);
}

template<class flt_t, class acc_t>
void ElectrodeVectorIntel::pair_contribution(IntelBuffers<flt_t,acc_t> *buffers,
  double *vector)
{
  const int inum = list->inum;
  const int nthreads = comm->nthreads;
  const int ago = neighbor->ago;
  const int nall = atom->nlocal + atom->nghost;

  if (buffers_stale) {
    if (_lrt == 0 && ago != 0 && fix->separate_buffers() == 0) {
      fix->start_watch(TIME_PACK);
      int packthreads;
      if (nthreads > INTEL_HTHREADS) packthreads = nthreads;
      else packthreads = 1;
      #if defined(_OPENMP)
      #pragma omp parallel shared(packthreads) if (packthreads > 1)
      #endif
      {
        int ifrom, ito, tid;
        IP_PRE_omp_range_id_align(ifrom, ito, tid, nall,
                                  packthreads, sizeof(ATOM_T));
        buffers->thr_pack(ifrom,ito,ago);
      }
      fix->stop_watch(TIME_PACK);
    }
    buffers_stale = false;
  }

  ATOM_T *_noalias const x = buffers->get_x(0);
  flt_t *_noalias const q = buffers->get_q(0);
  int nlocal = atom->nlocal;
  int nthr;
  if (nthreads > INTEL_HTHREADS) nthr = nthreads;
  else nthr = 1;
  int *mask = atom->mask;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;  // NOLINT
  const int newton_pair = force->newton_pair;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    const bool i_in_sensor = (mask[i] & groupbit);
    const bool i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (!(i_in_sensor || i_in_source)) continue;
    flt_t const xtmp = x[i].x;
    flt_t const ytmp = x[i].y;
    flt_t const ztmp = x[i].z;
    double const qi = (double) q[i];
    double const eta_i = etaflag ? atom->dvector[eta_index][i] : eta;
    int const itype = x[i].w;
    int const * _noalias const jlist = firstneigh[i];
    int jnum = numneigh[i];
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma ivdep
    #endif
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      const bool j_in_sensor = (mask[j] & groupbit);
      const bool j_in_source = !!(mask[j] & source_grpbit) != invert_source;
      const bool compute_ij = i_in_sensor && j_in_source;
      const bool compute_ji = (newton_pair || j < nlocal) && (j_in_sensor && i_in_source);
      if (!(compute_ij || compute_ji)) continue;
      const flt_t delx = xtmp - x[j].x;
      const flt_t dely = ytmp - x[j].y;
      const flt_t delz = ztmp - x[j].z;
      const int jtype = IP_PRE_dword_index(x[j].w);
      const flt_t rsq = delx * delx + dely * dely + delz * delz;
      if (rsq >= cutsq[itype][jtype]) continue;
      const double eta_j = etaflag ? atom->dvector[eta_index][j] : eta;
      double etaij;
      if (i_in_sensor && j_in_sensor) {
        etaij = eta_i * eta_j / sqrt(eta_i * eta_i + eta_j * eta_j);
      } else if (i_in_sensor) {
        etaij = eta_i;
      } else {
        assert(j_in_sensor);
        etaij = eta_j;
      }
      const double r = (double) sqrt(rsq);
      const double rinv = 1.0 / r;
      double aij = rinv;
      aij *= ElectrodeMath::safe_erfc(g_ewald * r);
      aij -= ElectrodeMath::safe_erfc(etaij * r) * rinv;
      const double qj = (double) q[j];
      if (i_in_sensor) { vector[i] += aij * qj; }
      if (j_in_sensor && (!invert_source || !i_in_sensor)) { vector[j] += aij * qi; }
    }
  }
}
