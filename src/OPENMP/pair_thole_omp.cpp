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

#include "pair_thole_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix_drude.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairTholeOMP::PairTholeOMP(LAMMPS *lmp) :
  PairThole(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairTholeOMP::compute(int eflag, int vflag)
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
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else                    eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else                    eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else                    eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairTholeOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;
  const int * _noalias const drudetype = fix_drude->drudetype;
  const tagint * _noalias const drudeid = fix_drude->drudeid;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int * const * const firstneigh      = list->firstneigh;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];

    // only on core-drude pair
    if (drudetype[type[i]] == NOPOL_TYPE) continue;

    const int di = domain->closest_image(i, atom->map(drudeid[i]));
    if (di < 0) error->one(FLERR, "Drude partner of atom {} not found", atom->tag[i]);
    double qi;
    if (drudetype[type[i]] == DRUDE_TYPE)
      qi = q[i];
    else
      qi = -q[di];

    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const int itype    = type[i];
    const double * _noalias const cutsqi   = cutsq[itype];
    const double * _noalias const ascreeni = ascreen[itype];
    const double * _noalias const scalei   = scale[itype];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // only on core-drude pair, but not into the same pair
      if (drudetype[type[j]] == NOPOL_TYPE || j == di) continue;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq  = delx*delx + dely*dely + delz*delz;
      const int jtype   = type[j];

      if (rsq < cutsqi[jtype]) {

        // get dq of the core via the drude charge
        // look up the drude partner only within the cutoff, since it may
        // not be present as a ghost atom for more distant neighbors

        double qj;
        if (drudetype[type[j]] == DRUDE_TYPE)
          qj = q[j];
        else {
          const int dj = domain->closest_image(j, atom->map(drudeid[j]));
          if (dj < 0) error->one(FLERR, "Drude partner of atom {} not found", atom->tag[j]);
          qj = -q[dj];
        }

        const double r2inv  = 1.0/rsq;
        const double rinv   = sqrt(r2inv);
        const double r      = sqrt(rsq);
        const double asr    = ascreeni[jtype] * r;
        const double exp_asr = exp(-asr);
        const double dcoul  = qqrd2e * qi * qj * scalei[jtype] * rinv;
        const double factor_f = 0.5*(2. + (exp_asr * (-2. - asr * (2. + asr)))) - factor_coul;
        double factor_e = 0.0;
        if (EFLAG) factor_e = 0.5*(2. - (exp_asr * (2. + asr))) - factor_coul;
        const double fpair  = factor_f * dcoul * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        double ecoul = 0.0;
        if (EFLAG) ecoul = factor_e * dcoul;

        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                                  0.0, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairTholeOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairThole::memory_usage();
  return bytes;
}
