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

#include "pair_mesomem_dipole_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMesomemDipoleOMP::PairMesomemDipoleOMP(LAMMPS *lmp) :
    PairMesomemDipole(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairMesomemDipoleOMP::compute(int eflag, int vflag)
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
        if (force->newton_pair)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair)
        eval<0, 0, 1>(ifrom, ito, thr);
      else
        eval<0, 0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairMesomemDipoleOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, rsq, evdwl, factor_lj, inv_mag;
  double del[3], ni[3], fi[3], ti[3], tj[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *const *const torque = thr->get_torque();
  double *const *const mu = atom->mu;
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_lj = force->special_lj;
  double fxtmp, fytmp, fztmp, t1tmp, t2tmp, t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = t1tmp = t2tmp = t3tmp = 0.0;

    // particles with a zero dipole moment only have the isotropic interaction

    const double *nip = nullptr;
    if (mu[i][3] > 0.0) {
      inv_mag = 1.0 / mu[i][3];
      ni[0] = mu[i][0] * inv_mag;
      ni[1] = mu[i][1] * inv_mag;
      ni[2] = mu[i][2] * inv_mag;
      nip = ni;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      del[0] = xtmp - x[j].x;
      del[1] = ytmp - x[j].y;
      del[2] = ztmp - x[j].z;
      rsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];

      if (rsq < cutsq[itype][jtype]) {
        evdwl = mesomem_analytic(itype, jtype, rsq, del, nip, mu[j], fi, ti, tj);

        if (factor_lj != 1.0) {
          evdwl *= factor_lj;
          MathExtra::scale3(factor_lj, fi);
          MathExtra::scale3(factor_lj, ti);
          MathExtra::scale3(factor_lj, tj);
        }

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        t1tmp += ti[0];
        t2tmp += ti[1];
        t3tmp += ti[2];

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fi[0];
          f[j].y -= fi[1];
          f[j].z -= fi[2];
          torque[j][0] += tj[0];
          torque[j][1] += tj[1];
          torque[j][2] += tj[2];
        }

        if (EVFLAG)
          ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fi[0], fi[1], fi[2], del[0],
                           del[1], del[2], thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    torque[i][0] += t1tmp;
    torque[i][1] += t2tmp;
    torque[i][2] += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairMesomemDipoleOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMesomemDipole::memory_usage();

  return bytes;
}
