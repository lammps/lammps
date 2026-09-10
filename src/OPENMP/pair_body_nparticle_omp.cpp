// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "pair_body_nparticle_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBodyNparticleOMP::PairBodyNparticleOMP(LAMMPS *lmp) :
  PairBodyNparticle(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairBodyNparticleOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow per-atom data structures if needed

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->create(dnum, nall, "pair:dnum");
    memory->create(dfirst, nall, "pair:dfirst");
  }
  for (int i = 0; i < nall; i++) dnum[i] = 0;
  ndiscrete = 0;

  // pre-initialize body sub-particle data for all body atoms
  // must be done serially before the parallel force computation

  int *body = atom->body;
  for (int i = 0; i < nall; i++) {
    if (body[i] >= 0) body2space(i);
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
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
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBodyNparticleOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  int ni, nj, npi, npj, ifirst, jfirst;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcelj;
  double xi[3], xj[3], fi[3], fj[3], ti[3], tj[3];
  double *dxi, *dxj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  const double * const * _noalias const x = atom->x;
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  auto * _noalias const torque = (dbl3_t *) thr->get_torque()[0];
  const int * _noalias const body = atom->body;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // body sub-particle data has already been initialized

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      // body/body interactions = NxM sub-particles

      evdwl = 0.0;
      if (body[i] >= 0 && body[j] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];
        npj = dnum[j];
        jfirst = dfirst[j];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni];

          for (nj = 0; nj < npj; nj++) {
            dxj = discrete[jfirst+nj];

            xi[0] = x[i][0] + dxi[0];
            xi[1] = x[i][1] + dxi[1];
            xi[2] = x[i][2] + dxi[2];
            xj[0] = x[j][0] + dxj[0];
            xj[1] = x[j][1] + dxj[1];
            xj[2] = x[j][2] + dxj[2];

            delx = xi[0] - xj[0];
            dely = xi[1] - xj[1];
            delz = xi[2] - xj[2];
            rsq = delx*delx + dely*dely + delz*delz;

            r2inv = 1.0 / rsq;
            r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
            fpair = forcelj * r2inv;

            if (EFLAG)
              evdwl += r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);

            fi[0] = delx * fpair;
            fi[1] = dely * fpair;
            fi[2] = delz * fpair;
            f[i].x += fi[0];
            f[i].y += fi[1];
            f[i].z += fi[2];
            ti[0] = dxi[1]*fi[2] - dxi[2]*fi[1];
            ti[1] = dxi[2]*fi[0] - dxi[0]*fi[2];
            ti[2] = dxi[0]*fi[1] - dxi[1]*fi[0];
            torque[i].x += ti[0];
            torque[i].y += ti[1];
            torque[i].z += ti[2];

            if (NEWTON_PAIR || j < nlocal) {
              fj[0] = -delx * fpair;
              fj[1] = -dely * fpair;
              fj[2] = -delz * fpair;
              f[j].x += fj[0];
              f[j].y += fj[1];
              f[j].z += fj[2];
              tj[0] = dxj[1]*fj[2] - dxj[2]*fj[1];
              tj[1] = dxj[2]*fj[0] - dxj[0]*fj[2];
              tj[2] = dxj[0]*fj[1] - dxj[1]*fj[0];
              torque[j].x += tj[0];
              torque[j].y += tj[1];
              torque[j].z += tj[2];
            }
          }
        }

      // body/particle interaction = Nx1 sub-particles

      } else if (body[i] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni];

          xi[0] = x[i][0] + dxi[0];
          xi[1] = x[i][1] + dxi[1];
          xi[2] = x[i][2] + dxi[2];
          xj[0] = x[j][0];
          xj[1] = x[j][1];
          xj[2] = x[j][2];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          r2inv = 1.0 / rsq;
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fpair = forcelj * r2inv;

          if (EFLAG)
            evdwl += r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);

          fi[0] = delx * fpair;
          fi[1] = dely * fpair;
          fi[2] = delz * fpair;
          f[i].x += fi[0];
          f[i].y += fi[1];
          f[i].z += fi[2];
          ti[0] = dxi[1]*fi[2] - dxi[2]*fi[1];
          ti[1] = dxi[2]*fi[0] - dxi[0]*fi[2];
          ti[2] = dxi[0]*fi[1] - dxi[1]*fi[0];
          torque[i].x += ti[0];
          torque[i].y += ti[1];
          torque[i].z += ti[2];

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] = -delx * fpair;
            fj[1] = -dely * fpair;
            fj[2] = -delz * fpair;
            f[j].x += fj[0];
            f[j].y += fj[1];
            f[j].z += fj[2];
          }
        }

      // particle/body interaction = 1xN sub-particles

      } else if (body[j] >= 0) {
        npj = dnum[j];
        jfirst = dfirst[j];

        for (nj = 0; nj < npj; nj++) {
          dxj = discrete[jfirst+nj];

          xi[0] = x[i][0];
          xi[1] = x[i][1];
          xi[2] = x[i][2];
          xj[0] = x[j][0] + dxj[0];
          xj[1] = x[j][1] + dxj[1];
          xj[2] = x[j][2] + dxj[2];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          delz = xi[2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          r2inv = 1.0 / rsq;
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fpair = forcelj * r2inv;

          if (EFLAG)
            evdwl += r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);

          fi[0] = delx * fpair;
          fi[1] = dely * fpair;
          fi[2] = delz * fpair;
          f[i].x += fi[0];
          f[i].y += fi[1];
          f[i].z += fi[2];

          if (NEWTON_PAIR || j < nlocal) {
            fj[0] = -delx * fpair;
            fj[1] = -dely * fpair;
            fj[2] = -delz * fpair;
            f[j].x += fj[0];
            f[j].y += fj[1];
            f[j].z += fj[2];
            tj[0] = dxj[1]*fj[2] - dxj[2]*fj[1];
            tj[1] = dxj[2]*fj[0] - dxj[0]*fj[2];
            tj[2] = dxj[0]*fj[1] - dxj[1]*fj[0];
            torque[j].x += tj[0];
            torque[j].y += tj[1];
            torque[j].z += tj[2];
          }
        }

      // particle/particle interaction = 1x1 particles

      } else {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = forcelj * r2inv;

        if (EFLAG)
          evdwl += r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);

        f[i].x += delx * fpair;
        f[i].y += dely * fpair;
        f[i].z += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }
      }

      if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                               evdwl, 0.0, fpair, delx, dely, delz, thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBodyNparticleOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBodyNparticle::memory_usage();
  return bytes;
}
