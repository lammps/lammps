/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_edip_multi_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include "fmt/format.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathExtra;

// max number of interaction per atom for f(Z) environment potential

static constexpr int leadDimInteractionList = 64;

static inline void costheta_d(const double *dr_ij, const double r_ij, const double *dr_ik,
                              const double r_ik, double *dri, double *drj, double *drk)
{
  const double costheta = dot3(dr_ij, dr_ik) / r_ij / r_ik;
  scaleadd3(1 / r_ij / r_ik, dr_ik, -costheta / r_ij / r_ij, dr_ij, drj);
  scaleadd3(1 / r_ij / r_ik, dr_ij, -costheta / r_ik / r_ik, dr_ik, drk);
  scaleadd3(-1, drj, -1, drk, dri);
}

/* ---------------------------------------------------------------------- */

PairEDIPMultiOMP::PairEDIPMultiOMP(LAMMPS *lmp) : PairEDIPMulti(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEDIPMultiOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        eval<1, 1>(ifrom, ito, thr);
      } else {
        eval<1, 0>(ifrom, ito, thr);
      }
    } else
      eval<0, 0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region

  error_thr();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairEDIPMultiOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, k, ii, jj, kk, jnum;
  int itype, jtype, ktype, ijparam, ikparam, ijkparam;
  double xtmp, ytmp, ztmp, evdwl;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int preForceCoord_counter;

  double zeta_i;
  double dzetair;
  double fpair;
  double costheta;
  double dpairZ, dtripleZ;

  evdwl = 0.0;

  const int tid = thr->get_tid();

  // per-thread slice of the scratch array for the coordination forces

  double *pre_thrForceCoord = preForceCoord + 5 * tid * leadDimInteractionList;

  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    zeta_i = 0.0;
    int numForceCoordPairs = 0;

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;

    // all the neighbors of atom i

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // the per-neighbor scratch array has a fixed size

    if (check_error_thr((jnum > leadDimInteractionList), tid, FLERR,
                        fmt::format("Pair style edip/multi supports at most {} neighbors per "
                                    "atom; atom {} has {}",
                                    leadDimInteractionList, atom->tag[i], jnum)))
      return;

    // pre-loop to compute environment coordination f(Z)

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      double delx, dely, delz, r_ij;

      delx = x[j].x - xtmp;
      dely = x[j].y - ytmp;
      delz = x[j].z - ztmp;
      r_ij = delx * delx + dely * dely + delz * delz;

      jtype = map[type[j]];
      const Param &param = params[elem3param[itype][jtype][jtype]];

      if (r_ij > param.cutsq) continue;

      r_ij = sqrt(r_ij);

      // zeta and its derivative dZ/dr

      if (r_ij < param.cutoffC)
        zeta_i += 1.0;
      else {
        double fc, fdr;
        edip_fc(r_ij, param, fc, fdr);
        zeta_i += fc;
        dzetair = -fdr / r_ij;

        preForceCoord_counter = numForceCoordPairs * 5;
        pre_thrForceCoord[preForceCoord_counter + 0] = dzetair;
        pre_thrForceCoord[preForceCoord_counter + 1] = delx;
        pre_thrForceCoord[preForceCoord_counter + 2] = dely;
        pre_thrForceCoord[preForceCoord_counter + 3] = delz;
        pre_thrForceCoord[preForceCoord_counter + 4] = j;
        numForceCoordPairs++;
      }
    }

    // two-body interactions

    dpairZ = 0;
    dtripleZ = 0;

    for (jj = 0; jj < jnum; jj++) {
      double dr_ij[3], r_ij, f_ij[3];

      j = jlist[jj];
      j &= NEIGHMASK;

      dr_ij[0] = x[j].x - xtmp;
      dr_ij[1] = x[j].y - ytmp;
      dr_ij[2] = x[j].z - ztmp;
      r_ij = dr_ij[0] * dr_ij[0] + dr_ij[1] * dr_ij[1] + dr_ij[2] * dr_ij[2];

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (r_ij > params[ijparam].cutsq) continue;

      r_ij = sqrt(r_ij);

      // potential energy and force
      // since pair i-j is different from pair j-i, double counting is
      // already considered in constructing the potential

      double fdr, fdZ;
      edip_pair(r_ij, zeta_i, params[ijparam], evdwl, fdr, fdZ);
      fpair = -fdr / r_ij;
      dpairZ += fdZ;

      f[i].x -= fpair * dr_ij[0];
      f[i].y -= fpair * dr_ij[1];
      f[i].z -= fpair * dr_ij[2];

      f[j].x += fpair * dr_ij[0];
      f[j].y += fpair * dr_ij[1];
      f[j].z += fpair * dr_ij[2];

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, newton_pair, evdwl, 0.0, fpair, -dr_ij[0], -dr_ij[1],
                     -dr_ij[2], thr);

      // three-body Forces

      for (kk = jj + 1; kk < jnum; kk++) {
        double dr_ik[3], r_ik, f_ik[3];

        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        dr_ik[0] = x[k].x - xtmp;
        dr_ik[1] = x[k].y - ytmp;
        dr_ik[2] = x[k].z - ztmp;
        r_ik = dr_ik[0] * dr_ik[0] + dr_ik[1] * dr_ik[1] + dr_ik[2] * dr_ik[2];

        if (r_ik > params[ikparam].cutsq) continue;

        r_ik = sqrt(r_ik);

        costheta = dot3(dr_ij, dr_ik) / r_ij / r_ik;

        double v1, v2, v3, v4, v5, v6, v7;

        edip_fcut3(r_ij, params[ijparam], v1, v2);
        edip_fcut3(r_ik, params[ikparam], v3, v4);
        edip_h(costheta, zeta_i, params[ijkparam], v5, v6, v7);

        // potential energy and forces
        evdwl = v1 * v3 * v5;
        dtripleZ += v1 * v3 * v7;

        double dri[3], drj[3], drk[3];
        double dhl, dfr;

        dhl = v1 * v3 * v6;

        costheta_d(dr_ij, r_ij, dr_ik, r_ik, dri, drj, drk);

        f_ij[0] = -dhl * drj[0];
        f_ij[1] = -dhl * drj[1];
        f_ij[2] = -dhl * drj[2];
        f_ik[0] = -dhl * drk[0];
        f_ik[1] = -dhl * drk[1];
        f_ik[2] = -dhl * drk[2];

        dfr = v2 * v3 * v5;
        fpair = -dfr / r_ij;

        f_ij[0] += fpair * dr_ij[0];
        f_ij[1] += fpair * dr_ij[1];
        f_ij[2] += fpair * dr_ij[2];

        dfr = v1 * v4 * v5;
        fpair = -dfr / r_ik;

        f_ik[0] += fpair * dr_ik[0];
        f_ik[1] += fpair * dr_ik[1];
        f_ik[2] += fpair * dr_ik[2];

        f[j].x += f_ij[0];
        f[j].y += f_ij[1];
        f[j].z += f_ij[2];

        f[k].x += f_ik[0];
        f[k].y += f_ik[1];
        f[k].z += f_ik[2];

        f[i].x -= f_ij[0] + f_ik[0];
        f[i].y -= f_ij[1] + f_ik[1];
        f[i].z -= f_ij[2] + f_ik[2];

        if (EVFLAG) ev_tally3_thr(this, i, j, k, evdwl, 0.0, f_ij, f_ik, dr_ij, dr_ik, thr);
      }
    }

    // forces due to environment coordination f(Z)
    for (int idx = 0; idx < numForceCoordPairs; idx++) {
      double delx, dely, delz;

      preForceCoord_counter = idx * 5;
      dzetair = pre_thrForceCoord[preForceCoord_counter + 0];
      delx = pre_thrForceCoord[preForceCoord_counter + 1];
      dely = pre_thrForceCoord[preForceCoord_counter + 2];
      delz = pre_thrForceCoord[preForceCoord_counter + 3];
      j = static_cast<int>(pre_thrForceCoord[preForceCoord_counter + 4]);

      dzetair *= (dpairZ + dtripleZ);

      f[j].x += dzetair * delx;
      f[j].y += dzetair * dely;
      f[j].z += dzetair * delz;

      f[i].x -= dzetair * delx;
      f[i].y -= dzetair * dely;
      f[i].z -= dzetair * delz;

      evdwl = 0.0;
      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, newton_pair, evdwl, 0.0, dzetair, -delx, -dely, -delz,
                     thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairEDIPMultiOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEDIPMulti::memory_usage();

  return bytes;
}
