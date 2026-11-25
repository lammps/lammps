// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#include "pair_lebedeva_z_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLebedevaZOMP::PairLebedevaZOMP(LAMMPS *lmp) :
  PairLebedevaZ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLebedevaZOMP::compute(int eflag, int vflag)
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
        if (vflag_either) eval<1,1,1>(ifrom, ito, thr);
        else              eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (vflag_either) eval<1,0,1>(ifrom, ito, thr);
        else              eval<1,0,0>(ifrom, ito, thr);
      }
    } else eval<0,0,0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int VFLAG_EITHER>
void PairLebedevaZOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int * const * const firstneigh      = list->firstneigh;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const int jtype = type[j];

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      // rho^2 = r^2 - z^2
      const double rhosq = delx*delx + dely*dely;
      const double rsq   = rhosq + delz*delz;

      if (rsq < cutsq[itype][jtype]) {

        const int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param& p = params[iparam_ij];

        const double r  = sqrt(rsq);
        const double r6 = rsq*rsq*rsq;
        const double r8 = r6*rsq;

        const double exp1 = exp(-p.alpha*(r-p.z0));
        const double exp2 = exp(-p.lambda1*rhosq);
        const double exp3 = exp(-p.lambda2*(delz*delz-p.z02));
        const double sumD = 1+p.D1*rhosq+p.D2*rhosq*rhosq;
        const double Ulm  = -p.A*p.z06/r6 + p.B*exp1 + p.C*sumD*exp2*exp3;

        // derivatives
        const double fpair = -6.0*p.A*p.z06/r8 + p.B*p.alpha*exp1/r;
        const double der   = p.D1 + 2*p.D2*rhosq - p.lambda1*sumD;
        const double fxy   = 2*p.C*exp2*exp3*der;
        const double fz    = 2*p.C*p.lambda2*sumD*exp2*exp3;

        f[i].x += delx*(fpair-fxy);
        f[i].y += dely*(fpair-fxy);
        f[i].z += delz*(fpair+fz);
        // newton_pair is always on (enforced by PairLebedevaZ::init_style)
        f[j].x -= delx*(fpair-fxy);
        f[j].y -= dely*(fpair-fxy);
        f[j].z -= delz*(fpair+fz);

        double evdwl = 0.0;
        if (EFLAG) evdwl = Ulm - offset[itype][jtype];

        if (EVFLAG) {
          ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1,
                       evdwl, 0.0, fpair, delx, dely, delz, thr);
          if (VFLAG_EITHER) {
            double fi[3], fj[3];
            fi[0] = -delx * fxy;
            fi[1] = -dely * fxy;
            fi[2] =  delz * fz;
            fj[0] =  delx * fxy;
            fj[1] =  dely * fxy;
            fj[2] = -delz * fz;
            v_tally2_newton_thr(this, i, fi, atom->x[i], thr);
            v_tally2_newton_thr(this, j, fj, atom->x[j], thr);
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLebedevaZOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLebedevaZ::memory_usage();
  return bytes;
}
