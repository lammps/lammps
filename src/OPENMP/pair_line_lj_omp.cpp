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

#include "pair_line_lj_omp.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLineLJOMP::PairLineLJOMP(LAMMPS *lmp) :
  PairLineLJ(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLineLJOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  int *line = atom->line;
  int *type = atom->type;

  // grow discrete list if necessary and initialize
  if (nall > nmax) {
    nmax = nall;
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->create(dnum,nall,"pair:dnum");
    memory->create(dfirst,nall,"pair:dfirst");
  }
  for (int i = 0; i < nall; i++) dnum[i] = 0;
  ndiscrete = 0;

  // pre-initialize all line atoms serially to avoid races on discrete[]
  // in the parallel region below
  for (int i = 0; i < nall; i++) {
    if (line[i] >= 0) discretize(i, subsize[type[i]]);
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
void PairLineLJOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,term1,term2,sig,sig3,forcelj;
  double xi[2],xj[2],fi[2],dxi,dxj,dyi,dyj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const double * _noalias const * _noalias const x = atom->x;
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  auto * _noalias const torque = (dbl3_t *) thr->get_torque()[0];
  const int * _noalias const line = atom->line;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

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
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      // line/line interactions = NxN particles

      evdwl = 0.0;
      fpair = 0.0;
      if (line[i] >= 0 && line[j] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];
        npj = dnum[j];
        jfirst = dfirst[j];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;

          for (nj = 0; nj < npj; nj++) {
            dxj = discrete[jfirst+nj].dx;
            dyj = discrete[jfirst+nj].dy;

            xi[0] = x[i][0] + dxi;
            xi[1] = x[i][1] + dyi;
            xj[0] = x[j][0] + dxj;
            xj[1] = x[j][1] + dyj;

            delx = xi[0] - xj[0];
            dely = xi[1] - xj[1];
            rsq = delx*delx + dely*dely;

            if (rsq >= cutsubsq[itype][jtype]) continue;

            sig = sigma[itype][jtype];
            sig3 = sig*sig*sig;
            term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
            term1 = 2.0 * term2 * sig3*sig3;
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (term1*r6inv - term2);
            fpair = forcelj*r2inv;

            if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

            fi[0] = delx*fpair;
            fi[1] = dely*fpair;

            f[i].x += fi[0];
            f[i].y += fi[1];
            torque[i].z += dxi*fi[1] - dyi*fi[0];

            if (NEWTON_PAIR || j < nlocal) {
              f[j].x -= fi[0];
              f[j].y -= fi[1];
              torque[j].z -= dxj*fi[1] - dyj*fi[0];
            }
          }
        }

      // line/particle interaction = Nx1 particles

      } else if (line[i] >= 0) {
        npi = dnum[i];
        ifirst = dfirst[i];

        for (ni = 0; ni < npi; ni++) {
          dxi = discrete[ifirst+ni].dx;
          dyi = discrete[ifirst+ni].dy;

          xi[0] = x[i][0] + dxi;
          xi[1] = x[i][1] + dyi;
          xj[0] = x[j][0];
          xj[1] = x[j][1];

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          rsq = delx*delx + dely*dely;

          if (rsq >= cutsubsq[itype][jtype]) continue;

          sig = sigma[itype][jtype];
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] = delx*fpair;
          fi[1] = dely*fpair;
          f[i].x += fi[0];
          f[i].y += fi[1];
          torque[i].z += dxi*fi[1] - dyi*fi[0];

          if (NEWTON_PAIR || j < nlocal) {
            f[j].x -= fi[0];
            f[j].y -= fi[1];
          }
        }

      // particle/line interaction = 1xN particles

      } else if (line[j] >= 0) {
        npj = dnum[j];
        jfirst = dfirst[j];

        for (nj = 0; nj < npj; nj++) {
          dxj = discrete[jfirst+nj].dx;
          dyj = discrete[jfirst+nj].dy;

          xi[0] = x[i][0];
          xi[1] = x[i][1];
          xj[0] = x[j][0] + dxj;
          xj[1] = x[j][1] + dyj;

          delx = xi[0] - xj[0];
          dely = xi[1] - xj[1];
          rsq = delx*delx + dely*dely;

          if (rsq >= cutsubsq[itype][jtype]) continue;

          sig = sigma[itype][jtype];
          sig3 = sig*sig*sig;
          term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
          term1 = 2.0 * term2 * sig3*sig3;
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (term1*r6inv - term2);
          fpair = forcelj*r2inv;

          if (EFLAG) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

          fi[0] = delx*fpair;
          fi[1] = dely*fpair;
          f[i].x += fi[0];
          f[i].y += fi[1];

          if (NEWTON_PAIR || j < nlocal) {
            f[j].x -= fi[0];
            f[j].y -= fi[1];
            torque[j].z -= dxj*fi[1] - dyj*fi[0];
          }
        }

      // particle/particle interaction = 1x1 particles

      } else {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = forcelj*r2inv;

        if (EFLAG)
          evdwl += r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

        f[i].x += delx*fpair;
        f[i].y += dely*fpair;
        f[i].z += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }
      }

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLineLJOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLineLJ::memory_usage();

  return bytes;
}
