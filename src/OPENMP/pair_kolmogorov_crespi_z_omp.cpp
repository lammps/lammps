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

#include "pair_kolmogorov_crespi_z_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairKolmogorovCrespiZOMP::PairKolmogorovCrespiZOMP(LAMMPS *lmp) :
  PairKolmogorovCrespiZ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairKolmogorovCrespiZOMP::compute(int eflag, int vflag)
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
void PairKolmogorovCrespiZOMP::eval(int iifrom, int iito, ThrData * const thr)
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
      // rho^2 = r^2 - (n,r) = r^2 - z^2
      const double rhosq = delx*delx + dely*dely;
      const double rsq   = rhosq + delz*delz;

      if (rsq < cutsq[itype][jtype]) {

        const int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param &p = params[iparam_ij];

        const double r     = sqrt(rsq);
        const double r6    = rsq*rsq*rsq;
        const double r8    = r6*rsq;
        const double rdsq  = rhosq * p.delta2inv;    // (rho/delta)^2

        const double exp1  = exp(-p.lambda*(r-p.z0));
        const double exp2  = exp(-rdsq);

        const double sumC   = p.C0 + p.C2*rdsq + p.C4*rdsq*rdsq;
        const double sumC2  = (2*p.C2 + 4*p.C4*rdsq)*p.delta2inv;
        const double frho   = exp2*sumC;
        const double sumCff = p.C + 2*frho;

        // derivatives
        const double fpair  = -6.0*p.A*p.z06/r8 + p.lambda*exp1/r*sumCff;
        const double fpair1 = exp1*exp2*(4.0*p.delta2inv*sumC - 2.0*sumC2);
        const double fsum   = fpair + fpair1;

        f[i].x += delx * fsum;
        f[i].y += dely * fsum;
        // fi_z does not contain contributions from df/dr
        // because rho_ij does not depend on z_i or z_j
        f[i].z += delz * fpair;
        // newton_pair is always on (enforced by PairKolmogorovCrespiZ::init_style)
        f[j].x -= delx * fsum;
        f[j].y -= dely * fsum;
        f[j].z -= delz * fpair;

        double evdwl = 0.0;
        if (EFLAG) evdwl = -p.A*p.z06/r6 + exp1*sumCff - offset[itype][jtype];

        if (EVFLAG) {
          ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1,
                       evdwl, 0.0, fpair, delx, dely, delz, thr);
          if (VFLAG_EITHER) {
            double fi[3], fj[3];
            fi[0] =  delx * fpair1;
            fi[1] =  dely * fpair1;
            fi[2] =  0;
            fj[0] = -delx * fpair1;
            fj[1] = -dely * fpair1;
            fj[2] =  0;
            v_tally2_newton_thr(this, i, fi, atom->x[i], thr);
            v_tally2_newton_thr(this, j, fj, atom->x[j], thr);
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairKolmogorovCrespiZOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairKolmogorovCrespiZ::memory_usage();
  return bytes;
}
