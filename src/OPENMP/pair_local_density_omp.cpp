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

#include "pair_local_density_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLocalDensityOMP::PairLocalDensityOMP(LAMMPS *lmp) :
  PairLocalDensity(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLocalDensityOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  const int newton_pair = force->newton_pair;

  // grow localrho and fp arrays if necessary.
  // localrho is allocated with nthreads*nmax to support per-thread
  // accumulation: thread tid writes to localrho[k][tid*nall + i]
  if (atom->nmax > nmax) {
    memory->destroy(localrho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(localrho, nLD, nthreads*nmax, "pairLD:localrho");
    memory->create(fp, nLD, nmax, "pairLD:fp");
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag,nall)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // zero per-thread localrho slice (nall elements per thread)
    for (int k = 0; k < nLD; k++)
      memset(&localrho[k][tid*nall], 0, nall * sizeof(double));

    // -------------------------------------------------------------------
    // Pass 1: compute per-thread local densities

    {
      int i,j,jj,jnum,itype,jtype,k;
      double xtmp,ytmp,ztmp,delx,dely,delz,rsq,phi;
      int *jlist;

      const double * _noalias const * _noalias const x = atom->x;
      const int * _noalias const type = atom->type;
      const int nlocal = atom->nlocal;
      int * const * const firstneigh = list->firstneigh;
      int * const numneigh = list->numneigh;
      int * const ilist_local = list->ilist;

      for (int ii = ifrom; ii < ito; ii++) {
        i = ilist_local[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          jtype = type[j];

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          for (k = 0; k < nLD; k++) {
            if (rsq < lowercutsq[k]) {
              phi = 1.0;
            } else if (rsq > uppercutsq[k]) {
              phi = 0.0;
            } else {
              phi = c0[k] + rsq*(c2[k] + rsq*(c4[k] + c6[k]*rsq));
            }
            localrho[k][tid*nall + i] += phi * b[k][jtype];
            if (newton_pair || j < nlocal)
              localrho[k][tid*nall + j] += phi * b[k][itype];
          }
        }
      }
    }

    sync_threads();

    thr->timer(Timer::PAIR);

    // reduce per-thread localrho contributions into localrho[k][0..nall-1]
    for (int k = 0; k < nLD; k++)
      data_reduce_thr(localrho[k], nall, nthreads, 1, tid);

    sync_threads();

    // reverse comm: master thread only
#if defined(_OPENMP)
#pragma omp master
#endif
    { if (newton_pair) comm->reverse_comm(this); }

    sync_threads();

    // -------------------------------------------------------------------
    // Pass 2: compute embedding function derivative fp and embedding energy

    {
      int i, k, m;
      double uLD, p, *coeff;
      const int * _noalias const type = atom->type;
      const int * _noalias const ilist_local = list->ilist;
      const int nlocal = atom->nlocal;

      for (int ii = ifrom; ii < ito; ii++) {
        i = ilist_local[ii];
        int itype = type[i];
        uLD = 0.0;

        for (k = 0; k < nLD; k++) {
          if (!a[k][itype]) continue;

          if (localrho[k][i] <= rho_min[k]) {
            coeff = frho_spline[k][0];
            fp[k][i] = coeff[2];
            uLD += a[k][itype] * (coeff[6] + fp[k][i]*(localrho[k][i] - rho_min[k]));
          } else if (localrho[k][i] >= rho_max[k]) {
            coeff = frho_spline[k][nrho-2];
            fp[k][i] = coeff[0] + coeff[1] + coeff[2];
            uLD += a[k][itype] * ((coeff[3]+coeff[4]+coeff[5]+coeff[6])
                                  + fp[k][i]*(localrho[k][i] - rho_max[k]));
          } else {
            p = (localrho[k][i] - rho_min[k]) / delta_rho[k];
            m = static_cast<int>(p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            fp[k][i] = (coeff[0]*p + coeff[1])*p + coeff[2];
            uLD += a[k][itype] * (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
          }
        }

        if (eflag)
          e_tally_thr(this, i, i, nlocal, /* newton_pair */ 1, uLD, 0.0, thr);
      }
    }

    sync_threads();

    // forward comm: master thread only
#if defined(_OPENMP)
#pragma omp master
#endif
    { comm->forward_comm(this); }

    sync_threads();

    // -------------------------------------------------------------------
    // Pass 3: compute forces

    if (evflag) {
      if (eflag) {
        if (newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else             eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else             eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else             eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLocalDensityOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,k,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double rsqinv,dphi,fpair;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const double * _noalias const * _noalias const x = atom->x;
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      fpair = 0.0;
      if (rsq < cutforcesq) {
        rsqinv = 1.0/rsq;
        for (k = 0; k < nLD; k++) {
          if (rsq >= lowercutsq[k] && rsq < uppercutsq[k]) {
            dphi = rsq * (2.0*c2[k] + rsq*(4.0*c4[k] + 6.0*c6[k]*rsq));
            fpair += -(a[k][itype]*b[k][jtype]*fp[k][i]
                      + a[k][jtype]*b[k][itype]*fp[k][j]) * dphi;
          }
        }
        fpair *= rsqinv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        // eng_vdwl already accumulated in pass 2; pass evdwl=0 here
        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                                 0.0, 0.0, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLocalDensityOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLocalDensity::memory_usage();
  // Extra memory for per-thread localrho
  bytes += (double)(comm->nthreads - 1) * nmax * nLD * sizeof(double);

  return bytes;
}
